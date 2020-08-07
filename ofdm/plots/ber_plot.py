import sys
sys.path.append('/usr/local/share/itpp')

import math
import pyitpp
import numpy as np

from scipy import special
from matplotlib import pyplot as plt

data = pyitpp.itload('ofdm.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-*')

data = pyitpp.itload('ofdm_tdl.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-x')

#ebn0 = np.power(10, np.divide(ebn0_dB, 10.0))
#plt.plot(ebn0_dB, 0.5*(1 - np.sqrt(np.divide(ebn0, (2+ebn0)))), '-o')

plt.title('OFDM - BER vs EbN0 plot')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.legend(['AWGN channel', 'TDL channel'])
plt.grid(True)
plt.show()
