from pymongo.mongo_client import MongoClient
from dotenv import load_dotenv
import os
from f_uniprot_metadata_retrieve import uniprot_id_to_info_dict
import subprocess
from Bio import SeqIO




load_dotenv()
MONGODB_URI = os.environ['MONGODB_URI']
client = MongoClient(MONGODB_URI)
ROAR_DB = client['ROAR-DB']
Rieske_Collection = ROAR_DB['Rieske']
Relative_Collection = ROAR_DB['RelativeGenes']


the_value = "N/A"
documents = Rieske_Collection.find({"PredictedCluster": the_value}, {"_id": 1})
for document in documents:


    if document:
        # Update the `_id` field
        new_id = document["_id"].split('|')[1]
        document["_id"] = new_id

        # Insert the updated document with the new `_id`
        try:
            Rieske_Collection.insert_one(document)
            print("Document inserted with new _id.")
            
            # Delete the original document
            Rieske_Collection.delete_one({"_id": original_id})
            print("Original document deleted.")
        except Exception as e:
            print(f"Error occurred: {e}")
    else:
        print("Document with the original _id not found.")



extra_info_dict = uniprot_id_to_info_dict(processed_result["Accession Name"].split('|')[1])