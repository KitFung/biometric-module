c++ library: Armadillo
`sudo apt install libarmadillo-dev`

# input:
- expression matrix
- count table
- meta data
  - how to clasify the sample to 2 group

# need preprocessing function. Gene differential
- tell the program how to classify the sample to 2 group
- use which columns and value as one group
- Idea 1:
  - input columns list and list of value combination
    - all the sample meet any value combination will be classify as Group A
    - otherwise group as Group B

#
