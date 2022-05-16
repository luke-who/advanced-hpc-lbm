import matplotlib.pyplot as plt

x = list(range(1,29))

y_close_1_core_baseline = {"128x128":6.03,"256x256":41.927,"1024x1024":215.96}
y_close = {"128x128":[6.03,3.31,2.46,1.962,1.706,1.487,1.323,1.181,1.129,1.027,0.9772,0.912894,0.867844,0.849577,
                      0.907616,0.905401,0.885406,0.879234,0.831362,0.830440,0.845484,0.808501,0.826072,0.835150,0.802825,0.808165,0.830641,0.820046],
           "256x256":[41.927,22.098,15.63,12.11,10.31,8.76,7.61566,6.676421,6.14,5.558199,5.155802,4.802288,4.427067,4.229699,
                      4.066787,3.844124,3.897875,3.682617,3.512783,3.463155,3.29924,3.164147,3.191518,3.121992,2.964108,2.907392,2.835672,2.755821],
           "1024x1024":[215.96,112.02,78.354368,61.957343,52.482,46.177,42.283781,38.602632,38.20372,37.708683,38.0,37.656625,37.88184,35.103446,
                        29.580384,27.041985,25.394372,25.078169,24.154593,25.796074,24.153987,18.193484,18.36825,17.287242,14.878130,14.648641,10.847470,11.766496]}
speedup_close = {"128x128_close": [ y_close_1_core_baseline["128x128"]/y for y in y_close["128x128"] ],
                 "256x256_close": [ y_close_1_core_baseline["256x256"]/y for y in y_close["256x256"] ],
                 "1024x1024_close": [ y_close_1_core_baseline["1024x1024"]/y for y in y_close["1024x1024"] ]}

y_spread_1_core_baseline = {"128x128":6.059148,"256x256":42.342895,"1024x1024":216}
y_spread = {"128x128":[6.059148,3.44,2.480161,1.99,1.776327,1.549124,1.418207,1.296325,1.233494,1.175785,1.112508,1.065924,1.028584,1.004624,
                       0.968046,0.934015,0.951573,0.927530,0.896194,0.87423,0.848821,0.846980,0.869497,0.819434,0.807903,0.801849,0.819681,0.814029],
            "256x256":[42.342895,21.766366,15.122645,11.374259,9.896398,8.43675,7.474643,6.622923,6.198916,5.790573,5.356234,4.968588,4.716766,4.438499,
                       4.186567,3.977337,3.869471,3.741065,3.531481,3.476204,3.317721,3.246466,3.370546,3.046161,2.987644,2.912377,2.851546,2.786261],
            "1024x1024":[216,98.972997,74.015203,53.066081,47.742511,36.987042,34.768546,29.478796,28.14531,25.108678,23.081051,20.461968,20.391231,19.436266,
                         18.313639,17.913499,17.181659,15.45543,16.045018,15.845361,17.562408,12.990233,16.448278,11.586166,15.036812,13.579105,13.625886,13.259012]}
speedup_spread = {"128x128_spread": [ y_spread_1_core_baseline["128x128"]/y for y in y_spread["128x128"] ],
                 "256x256_spread": [ y_spread_1_core_baseline["256x256"]/y for y in y_spread["256x256"] ],
                 "1024x1024_spread": [ y_spread_1_core_baseline["1024x1024"]/y for y in y_spread["1024x1024"] ]}

color=''

for k,v in speedup_close.items():
    if k=="128x128_close":
        color='c'
    elif k=="256x256_close":
        color='g'
    else:
        color='m'
    plt.plot(x,v,'o-',markersize=4,linewidth=1,color=color,label=k)
plt.legend(loc="upper left")
# plt.title('OMP_PROC_BIND=close')
plt.xlabel('#cores',fontsize=12)
plt.ylabel('speed-up',fontsize=12)
plt.tight_layout()
# plt.show() #this cannot be used together with plt.savefig() to produce the correct svg files
plt.savefig("svg/close.svg")
plt.close()

for k,v in speedup_spread.items():
    if k=="128x128_spread":
        color='c'
    elif k=="256x256_spread":
        color='g'
    else:
        color='m'
    plt.plot(x,v,marker='v',markersize=4,linestyle='-',linewidth=1,color=color,label=k)
plt.legend(loc="upper left")
# plt.title('OMP_PROC_BIND=spread')
plt.xlabel('#cores',fontsize=12)
plt.ylabel('speed-up',fontsize=12)
plt.tight_layout()
# plt.show() #this cannot be used together with plt.savefig() to produce the correct svg files
plt.savefig("svg/spread.svg")
plt.close()

