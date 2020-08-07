import sys
sys.path.append('/usr/local/share/itpp')

import math
import pyitpp
import numpy as np

from scipy import special
from matplotlib import pyplot as plt

fig, ax = plt.subplots()

data = pyitpp.itload('unquantized_feedback_rate_1by8.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-*')

data = pyitpp.itload('unquantized_feedback_rate_1by16.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-x')

data = pyitpp.itload('unquantized_feedback_rate_1by32.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-o')

textstr = "frame feedback rate - 1/8"
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.4, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.title('2x2 MIMO-OFDM Precoding (unquantized feedback)')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.legend(['subcarrier feedback rate - 1/8', 'subcarrier feedback rate - 1/16', 'subcarrier feedback rate - 1/32'], loc='lower left')
plt.grid(True)
plt.show()
#-------------------------------------------------------------------------------
fig, ax = plt.subplots()

data = pyitpp.itload('vector_quantization_2bit.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-*')

data = pyitpp.itload('vector_quantization_3bit.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-x')

data = pyitpp.itload('vector_quantization_4bit.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-o')

textstr = "frame feedback rate - 1/8\nsubcarrier feedback rate - 1/8"
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.4, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.title('2x2 MIMO-OFDM Precoding (vector quantization)')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.legend(['2bit quantizer', '3bit quantizer', '4bit quantizer'], loc='lower left')
plt.grid(True)
plt.show()

#-------------------------------------------------------------------------------
fig, ax = plt.subplots()

data = pyitpp.itload('unquantized_feedback_rate_1by8.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-*')

data = pyitpp.itload('vector_quantization_4bit.it')
ber = data.get('ber')
ebn0_dB = data.get('ebno_dB')
plt.plot(ebn0_dB, ber, '-o')

textstr = "frame feedback rate - 1/8\nsubcarrier feedback rate - 1/8"
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.4, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.title('2x2 MIMO-OFDM Precoding')
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('BER')
plt.yscale('log')
plt.legend(['unquantized', 'vector quantizer - 4bit'], loc='lower left')
plt.grid(True)
plt.show()
