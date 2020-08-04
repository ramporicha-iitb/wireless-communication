import sys
sys.path.append('/usr/local/share/itpp')

import math
import pyitpp
import numpy as np

from scipy import special
from matplotlib import pyplot as plt

data = pyitpp.itload('mimo_ofdm_tdl.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-*')

data = pyitpp.itload('alamouti/2x2_alamouti.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-x')

plt.title('2x2 MIMO-OFDM - BER vs EbN0 plot')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.legend(['Zero Forcing', 'Alamouti Coding'])
plt.grid(True)
plt.show()
