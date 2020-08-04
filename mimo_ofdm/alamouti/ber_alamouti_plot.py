import sys
sys.path.append('/usr/local/share/itpp')

import math
import pyitpp
import numpy as np

from scipy import special
from matplotlib import pyplot as plt

data = pyitpp.itload('2x1_alamouti.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-*')

data = pyitpp.itload('2x2_alamouti.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-x')

data = pyitpp.itload('2x3_alamouti.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-o')

plt.title('Alamouti Coding - BER vs EbN0 plot')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.legend(['2x1 MIMO', '2x2 MIMO', '2x3 MIMO'])
plt.grid(True)
plt.show()
