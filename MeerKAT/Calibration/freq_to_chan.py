"""
Author: Joppe Swart
Description: Simple file to get channels from frequency
"""

channels = 3723

freq_0 = 890.064
frequency = 1000

while frequency is not 'stop':
    print('From what frequency do you want the channel?')
    frequency = input()
    for i in range(channels):
       if frequency > freq_0 + i*0.208984 and frequency < freq_0 + (i+1)*0.208984:
           print('this corresponds to channel ', i)


