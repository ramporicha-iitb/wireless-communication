import sys
sys.path.append('/usr/local/share/itpp')

import math
import pyitpp
import numpy as np

from scipy import special
from matplotlib import pyplot as plt

fig, ax = plt.subplots()

data = pyitpp.itload('mimo_ofdm_tdl.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-*')

data = pyitpp.itload('alamouti/2x2_alamouti.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-^')

data = pyitpp.itload('precoding/unquantized_feedback_rate_1by8.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-o')

data = pyitpp.itload('precoding/vector_quantization_4bit.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-x')

textstr = "frame feedback rate - 1/8\nsubcarrier feedback rate - 1/8\nvector quantizer - 4bit"
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.01, 0.4, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.title('2x2 MIMO-OFDM - BER vs EbN0 plot')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.legend(['Zero Forcing', 'Alamouti Coding', 'Precoding (unquantized)', 'Precoding (vector quantization)'])
plt.grid(True)
plt.show()
