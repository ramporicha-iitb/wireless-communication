import sys
sys.path.append('/usr/local/share/itpp')

import math
import pyitpp
import numpy as np

from scipy import special
from matplotlib import pyplot as plt

data = pyitpp.itload('ofdm.it')
ber = data.get('ber')
ebno_dB = data.get('ebno_dB')
plt.plot(ebno_dB, ber)

plt.title('OFDM - BER vs EbN0 plot')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.show()
