from neo import io
from matplotlib import pyplot as plt
plt.ion()
r = io.RawBinarySignalIO( filename = '/data/data/130904_peter/amplipex/PeterP-2013-09-04_001.dat')
seg = r.read_segment(lazy = False, cascade = True,)

plt.figure()
plt.plot(seg.analogsignals[0])
print seg.analogsignals  