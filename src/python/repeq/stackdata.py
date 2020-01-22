import glob
import obspy
import matplotlib.pyplot as plt

EQdates=glob.glob('./HawaiiafterEQ_W/20*')

sumD=0
for EQdate in EQdates:
    sac=glob.glob(EQdate+'/waveforms/*RSDD*mseed')
    if sac==[]:
        continue
    print(sac)
    #read data
    D=obspy.read(sac[0])
    D[0].detrend('linear')
    D[0].normalize()
    sumD=sumD+D[0].data[:10000]
    t=D[0].times()[:10000]
    D.clear()

plt.plot(sumD,'k')
plt.show()
