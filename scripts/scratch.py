import yaml
import json  

# load a yaml file into dict()
stream = open('config.yaml', 'r')
config = yaml.load(stream, Loader=yaml.SafeLoader)

# load a json file into dict()
stream = open('config.json', 'r')
config = json.load(stream)

REF_VERSION = config["references"][REF_GENOME]["version"]

wildcards = dict()
wildcards={"assayType" : "ChIP-Seq",
            'reference_version' : REF_VERSION,
          'project' : "LR1807201",
          'runID' : "N08851_SK_LR1807201_SEQ",
          'cycle' : "G1"}