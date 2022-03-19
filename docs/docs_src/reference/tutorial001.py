import requests
import json
from sadie.reference import G3Error

url = "https://g3.jordanrwillis.com/api/v1/genes?source=imgt&segment=V&common=human&limit=5"
response = requests.get(url)
response_json = response.json()
if response.status_code != 200:
    raise G3Error("Error: " + str(response.status_code))
print(json.dumps(response_json, indent=4))
json.dump(response_json, open("human_v.json", "w"), indent=4)
