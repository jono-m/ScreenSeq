import nd2reader
import numpy as np
import matplotlib.pyplot as plt

intensities = []
with nd2reader.ND2Reader(r"C:\Users\jonoj\Documents\dropletdata\ControlPWM.nd2") as nd2file:
    for frame in nd2file.metadata['frames']:
        intensity = np.mean(nd2file.get_frame(frame))
        intensities.append(intensity)
secPerFrame = 30 / len(intensities)
startTrim = 50
endTrim = 30
intensities = intensities[startTrim:-endTrim]
intensities = 80 * (intensities - np.min(intensities)) / (np.max(intensities) - np.min(intensities)) + 10
dutyCycle = [(0, 10),
             (22, 20),
             (49, 30),
             (77, 40),
             (114, 50),
             (142, 60),
             (168, 70),
             (198, 80),
             (230, 90),
             (len(intensities), 90)]
plt.plot([x[0] for x in dutyCycle], [x[1] for x in dutyCycle], drawstyle="steps-post", color="black", linestyle="--")
plt.plot(intensities, color="orange")
plt.legend([r"Input Duty Cycle (%)", "Output Intensity (au)"])
plt.xlabel("Time (sec)")
plt.xticks([x[0] for x in dutyCycle], label=[x[0]*secPerFrame for x in dutyCycle])
plt.show()
