import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#reading and formating the dataframe
df = pd.read_csv("PM0.csv")
df = df.drop(df.index[[0, 1]])
df = df.values.tolist()

#here taking each column and creating their own numpy array for analysis
time, channel_a, channel_d = [], [], []
for i in range(len(df)):
    time.append(float(df[i][0]))
    channel_a.append(float(df[i][1]))
    channel_d.append(float(df[i][2]))
#time, channel_a, channel_d = np.array(time), np.array(channel_a), np.array(channel_d)



def time_duration(time, channel):
    threshold = 0
    temparray = []
    first_index = channel.index(max(channel))
    last_index = 0
    for i in range(0, len(time)):
       if channel[i] > threshold:                   
            last_index = channel.index(channel[i])

       #if still having trouble with +/- voltages, try this, it should catch the last index positive or negative before return to BG
       #if abs(channel[i]) > threshold: last_index = channel.index(channel[i])     
            
    print("This is the first and last index: "+str(first_index)+' '+str(last_index))
    print("this is the last index postion in time: "+str(time[last_index]))
    print("this is the first index postion in time: "+str(time[first_index]))


    temparray = np.array(temparray)
    return abs(first_index - last_index)





time_of_pulse = time_duration(time, channel_d)   #this returns a difference in indices, which is proportional to time

#print("This is the time of the pulse from channel_a: "+str(time_of_pulse))
#time_of_pulse = time_duration(time, channel_d)
#print("this is the time of the pulse from channel_d: "+str(time_of_pulse))



time, channel_a, channel_d = np.array(time), np.array(channel_a), np.array(channel_d)   #redundant?
fig, ax = plt.subplots()
ax.plot(time, channel_a, label="Channel_A")
ax.plot(time, channel_d*100, label="Channel_D")
plt.legend()

ay = ax.twinx()
ax.grid()
ax.set_xticks(np.linspace(min(time), max(time), 23))
ax.set_yticks(np.linspace(min(channel_a), max(channel_a), 25))

ay.set_yticks(np.linspace(min(channel_d), max(channel_d), 25))
ax.set_xlabel("Time (ms)")

ax.set_ylabel('mV: Channel_A')
ay.set_ylabel('V: Channel_D')



plt.show()

#print("This is the amplitude of channel_a: ", max(channel_a))
#print('This is the amplitude of channel_d: ', max(channel_d))
