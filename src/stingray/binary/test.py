from astropy.time import Time, TimeDelta

T_pi2 = 2455073.68504 - 2400000.5
Tnod = T_pi2 + Porb/2 

t = Time(58000, format='mjd')
print(t)
print(t.tdb)
print(t.unix)
print(t.value*86400)