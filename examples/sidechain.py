import glob
import multiprocessing
import pickle
from functools import reduce
from cpeptools import metrics
import mdtraj as md

name  = "W"
paths = [i.strip() for i in open("/fileserver/elm4/shuwang/water_models/{}/traj_order.txt".format(name))]
out_path = "/home/shuwang/Documents/Modelling/CP/Codes/"

print((paths))
cores = multiprocessing.cpu_count()

def wrapper(path) :return {path : metrics.identify_interference(md.load(path))}

with multiprocessing.Pool(processes=cores) as pool:
    out = pool.map(wrapper, paths )
out = reduce(lambda a, b : {**a , **b}, out)
pickle.dump(out, open(out_path + "{}.pickle".format(name), "wb"))
# print(out)


def wrapper(path) :return {path : list(metrics.calculate_eccentricity(md.load(path)))}

with multiprocessing.Pool(processes=cores) as pool:
    out = pool.map(wrapper, paths )
out = reduce(lambda a, b : {**a , **b}, out)
# print(out)
pickle.dump(out, open(out_path + "{}_eccen.pickle".format(name), "wb"))
