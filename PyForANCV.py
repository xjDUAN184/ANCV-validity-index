import matlab.engine
eng = matlab.engine.start_matlab()

"""
Python env: 3.10.9; matlab env: R2021a
Configuration method:
    1. Download the matlab client and keep your installation path
    2. Find the python engine in the matlab installation path: matlab\extern\engines\python
    3. Execute python setup.py install(cmd or Anaconda Prompt) under the path of python engine
Or you can look for tutorials on calling matlab programs through Python on the matlab official website
"""
K = 10 # k-nearest neighbors
SNN_thr = 3 # shared nearest neighbors
cluster_name_list = ["CTCEHC", "NTHC", "kmeans"]
cluster_name = cluster_name_list[0]
"""
If you are passing .txt file, you need to pass one argument: data_path
data_path is data
"""
# data_path = "./Aggregation.txt" # To modify the path here to your absolute path

"""
If you are passing .mat file, you need to pass in two arguments: data_path and data_group_path
data_path is data, data_group_path is label
"""
data_path = "./dataset/D12.mat" # To modify the path here to your absolute path
data_group_path = "./dataset/D12Group.mat" # To modify the path here to your absolute path

# Call a MATLAB program from Python script
eng.ANCV_index_py(cluster_name, K, SNN_thr, data_path, data_group_path)
