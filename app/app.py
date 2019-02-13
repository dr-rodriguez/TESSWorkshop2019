from flask import Flask, render_template, request, redirect, make_response
import os
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.models.widgets import Panel, Tabs
from bokeh.io import output_file, reset_output
import numpy as np
import pandas as pd
from bokeh.palettes import Spectral5 as color_palette
from astroquery.mast import Observations
import astropy.units as u
from astropy.coordinates import SkyCoord
from bokeh.layouts import column, widgetbox
import pymongo

app_portal = Flask(__name__)

# Connect to pymongo
db_name = os.environ.get('EXOMAST_NAME')
client = pymongo.MongoClient(os.environ.get('EXOMAST_MONGO'))
db = client[db_name]  # database
planets = db.planets  # collection


def parse_s_region(s_region):
    ra = []
    dec = []
    counter = 0
    for elem in s_region.strip().split():
        try:
            value = float(elem)
        except ValueError:
            continue
        if counter % 2 == 0:
            ra.append(value)
        else:
            dec.append(value)
        counter += 1

    return {'ra': ra, 'dec': dec}


def projection(lon, lat, use='hammer'):
    """
    Convert x,y to Aitoff or Hammer projection. Lat and Lon should be in radians. RA=lon, Dec=lat
    """
    # TODO: Figure out why Aitoff is failing

    # Note that np.sinc is normalized (hence the division by pi)
    if use.lower() == 'hammer':  # Hammer
        x = 2.0 ** 1.5 * np.cos(lat) * np.sin(lon / 2.0) / np.sqrt(1.0 + np.cos(lat) * np.cos(lon / 2.0))
        y = np.sqrt(2.0) * np.sin(lat) / np.sqrt(1.0 + np.cos(lat) * np.cos(lon / 2.0))
    else:  # Aitoff, not yet working
        alpha_c = np.arccos(np.cos(lat) * np.cos(lon / 2.0))
        x = 2.0 * np.cos(lat) * np.sin(lon) / np.sinc(alpha_c / np.pi)
        y = np.sin(lat) / np.sinc(alpha_c / np.pi)
    return x, y


def make_sky_plot(proj='hammer'):
    """
    Make a sky plot and return a Bokeh figure
    Adapted from: https://github.com/astrocatalogs/astrocats/blob/master/scripts/hammertime.py#L93-L132 and the Open Supernova Catalog (https://sne.space/statistics/sky-locations/)
    """
    # source = ColumnDataSource(data=data)

    tools = "save,pan,wheel_zoom,box_zoom,reset"
    p = figure(tools=tools, title='', plot_width=800, plot_height=600,
               x_range=(-1.05 * (2.0 ** 1.5), 1.3 * 2.0 ** 1.5), y_range=(-2.0 * np.sqrt(2.0), 1.2 * np.sqrt(2.0)),
               min_border=0, min_border_bottom=0)

    # Initial figure formatting
    p.axis.visible = None
    p.outline_line_color = None
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None

    # Add the grid
    pi = np.pi
    rangepts = 50
    raseps = 12
    decseps = 12
    rarange = [-pi + i * 2.0 * pi / rangepts for i in range(0, rangepts + 1)]
    decrange = [-pi / 2.0 + i * pi / rangepts for i in range(0, rangepts + 1)]
    ragrid = [-pi + i * 2.0 * pi / raseps for i in range(0, raseps + 1)]
    decgrid = [-pi / 2.0 + i * pi / decseps for i in range(0, decseps + 1)]

    raxs = []
    rays = []
    for rg in ragrid:
        t1 = [projection(rg, x, use=proj) for x in decrange]
        tx, ty = zip(*t1)
        raxs.append(tx)
        rays.append(ty)

    decxs = []
    decys = []
    for dg in decgrid:
        t1 = [projection(x, dg, use=proj) for x in rarange]
        tx, ty = zip(*t1)
        decxs.append(tx)
        decys.append(ty)

    p.multi_line(raxs, rays, color='#bbbbbb')
    p.multi_line(decxs, decys, color='#bbbbbb')

    # Add the data
    # p.scatter('x', 'y', source=source, size=8, alpha=0.6)
    # tooltip = [("obs_id", "@obs_id")]
    # p.add_tools(HoverTool(tooltips=tooltip))

    # When clicked, go to the Summary page
    # url = "summary/@id"
    # taptool = p.select(type=TapTool)
    # taptool.callback = OpenURL(url=url)

    return p


def add_patches(p, obsDF, proj='hammer', maptype='equatorial', tooltip=[("obs_id", "@obs_id")]):
    for sec, color in zip([1, 2, 3, 4, 5], color_palette):
        ind = obsDF['sequence_number'] == sec

        # Add patches with the observation footprints
        patch_xs = [c['ra'] for c in obsDF['coords'][ind]]
        patch_ys = [c['dec'] for c in obsDF['coords'][ind]]

        data = {'x': patch_xs, 'y': patch_ys, 'obs_id': obsDF['obs_id'][ind]}

        # Project coordinates
        new_x, new_y, new_l, new_b = [], [], [], []
        for i in range(len(patch_ys)):
            c = SkyCoord(ra=patch_xs[i] * u.degree, dec=patch_ys[i] * u.degree)
            if maptype == 'ecliptic':
                temp_x, temp_y = projection(c.geocentrictrueecliptic.lon.radian - np.pi,
                                            c.geocentrictrueecliptic.lat.radian, use=proj)
            elif maptype == 'galactic':
                temp_x, temp_y = projection(c.galactic.l.radian - np.pi, c.galactic.b.radian, use=proj)
            else:  # default is equatorial
                temp_x, temp_y = projection(c.ra.radian - np.pi, c.dec.radian, use=proj)

            # Fix for cases that cross the meridian
            sign = 1 if np.mean(temp_x) >= 0 else -1
            temp_x = [abs(x) * sign for x in temp_x]

            new_x.append(temp_x)
            new_y.append(temp_y)
            # print(sec, temp_x, temp_y)

        data['x'], data['y'] = new_x, new_y

        p.patches('x', 'y', source=data, legend='Sector {}'.format(sec),
                  fill_color=color, fill_alpha=0.2, line_color="black", line_width=0.5)

    p.legend.click_policy = "hide"

    # Add hover tooltip for MAST observations
    if tooltip:
        p.add_tools(HoverTool(tooltips=tooltip))

    return p


def add_points(p, data, proj='hammer', maptype='equatorial',
               tooltip=[('Planet Name', '@planet_name'),("(RA, Dec)", "(@ra, @dec)"),
                        ('Catalog', '@catalog_name')]):

    if proj:
        c = SkyCoord(ra=np.array(data['ra']) * u.degree, dec=np.array(data['dec']) * u.degree)
        if maptype == 'ecliptic':
            data['x'], data['y'] = projection(c.geocentrictrueecliptic.lon.radian - np.pi,
                                        c.geocentrictrueecliptic.lat.radian, use=proj)
        elif maptype == 'galactic':
            data['x'], data['y'] = projection(c.galactic.l.radian - np.pi, c.galactic.b.radian, use=proj)
        else:
            data['x'], data['y'] = projection(c.ra.radian - np.pi, c.dec.radian, use=proj)

    source = ColumnDataSource(data=data)
    p.scatter('x', 'y', source=source, size=4, alpha=0.2, legend='ExoMast Sources')
    if tooltip:
        p.add_tools(HoverTool(tooltips=tooltip))

    p.legend.click_policy = "hide"


# Redirect to the main page
@app_portal.route('/')
@app_portal.route('/index')
@app_portal.route('/index.html')
def app_home():
    return render_template('index.html')


@app_portal.route('/tessffi', methods=['GET', 'POST'])
def app_tessffi():
    obsTable = Observations.query_criteria(dataproduct_type=["image"], obs_collection='TESS')
    obsDF = obsTable.to_pandas()
    obsDF['coords'] = obsDF.apply(lambda x: parse_s_region(x['s_region']), axis=1)

    p1 = make_sky_plot()
    p1 = add_patches(p1, obsDF, maptype='equatorial')
    p2 = make_sky_plot()
    p2 = add_patches(p2, obsDF, maptype='galactic')
    p3 = make_sky_plot()
    p3 = add_patches(p3, obsDF, maptype='ecliptic')

    tab1 = Panel(child=p1, title="Equatorial")
    tab2 = Panel(child=p2, title="Galatic")
    tab3 = Panel(child=p3, title="Ecliptic")
    tabs = Tabs(tabs=[tab1, tab2, tab3])

    script, div = components(tabs)

    return render_template('tessffi.html', script=script, plot=div)


@app_portal.route('/tessexomast', methods=['GET', 'POST'])
def app_tessexomast():
    global planets

    # Load data
    cursor = planets.aggregate([{'$group': {
        '_id': '$exoplanet_id',
        'planet_name': {'$push': '$planet_name'},
        'catalog_name': {'$push': '$catalog_name'},
        'ra': {'$push': '$ra'},
        'dec': {'$push': '$dec'}
    }},
        {'$project': {'_id': 1, 'planet_name': 1, 'catalog_name': 1,
                      'ra': {'$slice': ['$ra', 0, 1]},
                      'dec': {'$slice': ['$dec', 0, 1]},
                      }}])
    data_all = list(cursor)
    # for i in range(len(data_all)): print(data_all[i])

    df = pd.DataFrame(data_all)
    df['ra'] = pd.to_numeric(df['ra'].apply(lambda x: x[0]))
    df['dec'] = pd.to_numeric(df['dec'].apply(lambda x: x[0]))

    # Prepare TESS FFI
    obsTable = Observations.query_criteria(dataproduct_type=["image"], obs_collection='TESS')
    obsDF = obsTable.to_pandas()
    obsDF['coords'] = obsDF.apply(lambda x: parse_s_region(x['s_region']), axis=1)

    p1 = make_sky_plot()
    # add_patches(p1, obsDF, maptype='equatorial', tooltip=None)
    add_points(p1, df, maptype='equatorial')
    # p2 = make_sky_plot()
    # add_patches(p2, obsDF, maptype='galactic', tooltip=None)
    # add_points(p2, df, maptype='galactic')
    # p3 = make_sky_plot()
    # add_patches(p3, obsDF, maptype='ecliptic', tooltip=None)
    # add_points(p3, df, maptype='ecliptic')

    # tab1 = Panel(child=p1, title="Equatorial")
    # tab2 = Panel(child=p2, title="Galatic")
    # tab3 = Panel(child=p3, title="Ecliptic")
    # tabs = Tabs(tabs=[tab1])

    script, div = components(p1)

    return render_template('tessffi.html', script=script, plot=div)


@app_portal.route('/exomast', methods=['GET', 'POST'])
def app_exomast():

    # Get nexsci data
    cursor = planets.find({'catalog_name': 'nexsci'},
                          {'_id': 0, 'orbital_period.value': 1, 'planet_radius.value': 1})
    df_nexsci = list(cursor)
    cursor = planets.find({'catalog_name': 'exoplanets.org'},
                          {'_id': 0, 'orbital_period.value': 1, 'planet_radius.value': 1})
    df_exoplanets = list(cursor)
    cursor = planets.find({'catalog_name': 'koi'},
                          {'_id': 0, 'orbital_period.value': 1, 'planet_radius.value': 1})
    df_koi = list(cursor)
    cursor = planets.find({'catalog_name': 'TESS-DV'},
                          {'_id': 0, 'orbital_period.value': 1, 'planet_radius.value': 1})
    df_tessdv = list(cursor)



