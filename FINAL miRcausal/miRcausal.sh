rm mirfluence.output coverage_inference.txt test.edgelist mapping.json
python mapper.py "$1"
miRfluence/InfluenceModels -c input.txt
python map_back.py

