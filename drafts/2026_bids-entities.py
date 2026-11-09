"""
Build the ordered list of BIDS filename entity labels (e.g. 'sub', 'ses', 'task', ...)
by combining:
  - src/schema/rules/entities.yaml   -> defines the ORDER (long/internal keys)
  - src/schema/objects/entities.yaml -> defines each entity's short "name" (the label
                                          actually used in filenames, e.g. sub, ses, acq)
 
Requires: pip install pyyaml requests
"""
 
import requests
import yaml
 
RULES_URL = (
    "https://raw.githubusercontent.com/bids-standard/"
    "bids-specification/master/src/schema/rules/entities.yaml"
)
OBJECTS_URL = (
    "https://raw.githubusercontent.com/bids-standard/"
    "bids-specification/master/src/schema/objects/entities.yaml"
)
 
 
def fetch_yaml(url: str):
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return yaml.safe_load(resp.text)
 
 
def build_bids_label_order():
    # 1. Order of entities, using their long/internal keys
    #    e.g. ['subject', 'template', 'session', 'cohort', 'sample', 'task', ...]
    order_keys = fetch_yaml(RULES_URL)
 
    # 2. Entity definitions: internal key -> {'name': short_label, ...}
    #    e.g. {'subject': {'name': 'sub', ...}, 'ceagent': {'name': 'ce', ...}, ...}
    entity_defs = fetch_yaml(OBJECTS_URL)
 
    # 3. Map each ordered internal key to its short filename label
    bids_label_order = [entity_defs[key]["name"] for key in order_keys]
 
    return bids_label_order
 
 
if __name__ == "__main__":
    bids_label_order = build_bids_label_order()
    print(bids_label_order)
 
