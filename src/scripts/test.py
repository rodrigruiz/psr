import pickle

with open('/home/wecapstor3/capn/mppi148h/gtis/Antares_054890_gtis.pkl', 'rb') as f:
    gtis_file = pickle.load(f)
print(gtis_file)