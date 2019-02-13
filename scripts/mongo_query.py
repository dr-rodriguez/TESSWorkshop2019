# Hack to load mongoDB information

import os
import pymongo
import pandas as pd

db_name = os.getenv('EXOMAST_NAME')
client = pymongo.MongoClient(os.getenv('EXOMAST_MONGO'))
db = client[db_name]  # database
planets = db.planets  # collection

# Count entries
planets.count_documents({})

# Group by exoplanet_id
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
for i in range(len(data_all)): print(data_all[i])

df = pd.DataFrame(data_all)
df['ra'] = df['ra'].apply(lambda x: x[0])
df['dec'] = df['dec'].apply(lambda x: x[0])

