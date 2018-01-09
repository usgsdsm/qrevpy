#-------------------------------------------------Begin Scratch Curve Fit Code -------------------------------
# import numpy as np
# from numpy.linalg import lstsq
# from scipy.optimize import curve_fit
# import matplotlib.pyplot as plt
# import uncertainties as unc
# import uncertainties.unumpy as unp
# 
# 
# 
# x = np.array([2,3,4,5,6])
# y = np.array([3,17,28,40,51])
# func = lambda x, a, b: a * x**b
# bounds = [[-np.inf,0.01],[np.inf, 1]]
# popt, pcov = curve_fit(func, x, y, p0=[2, 1./6], bounds = bounds)
# 
# 
#  
# from scipy.stats import t
# 
# t_val1 = t.cdf(.025, 4)
# # t_val2 = t.ppf(.025, 4)
# 
# alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
# 
# n = len(y)    # number of data points
# p = len(popt) # number of parameters
# 
# ss_tot = np.sum((y - np.mean(y))**2)
# ss_res = np.sum((y - func(x, *popt))**2)
# r_squared = 1 - (ss_res/ss_tot)
# 
# upper, lower = [], []
# for j in range(len(popt)):
#     if bounds[0][j] == -np.inf and bounds[1][j] == np.inf:
#         lower.append(popt[j] - t_val1 * np.sqrt(np.diag(pcov)[j])) 
#         upper.append(popt[j] + t_val1 * np.sqrt(np.diag(pcov)[j])) 
#     else:
#         lower.append(popt[j])
#         upper.append(popt[j])
#     
# print(upper, lower)
# plt.plot(x, func(x, *popt), alpha=.5)
#  
# plt.plot(x, func(x, *upper))
# 
# plt.plot(x, func(x, *lower))
# plt.show()
#-----------------------------------------------------------End Scratch Curve Fit Code


# m, c = lstsq(np.array([x, [1,1,1,1,1]]).T, y_0)[0]
# print(m, c)
# import matplotlib.pyplot as plt
# plt.plot(x, y, 'o', label='Original data', markersize=10)
# plt.plot(x, m*x + c, 'r', label='Fitted line')
# plt.legend()
# plt.show()

# #DOF n - p in this case because basis functions are independent
# DOF = 4
# reduced_chi_squared_min = np.sum((y - func(x, *popt)) **2) / DOF
# 
# def confidence_band(x, dfdp, alpha, model, abswei = True):
#    
#     
#     C = pcov
#     n = 2
#     p = popt
#     N = len(x)
#     
#     if abswei:
#         cov_scale = 1.0
#     else:
#         cov_scale = reduced_chi_squared_min
#         
#     df2 = np.zeros(N)
#     for j in range(n):
#         for k in range(n):
#             df2 += dfdp[j](x , *p) * dfdp[k](x, *p) * C[j,k]
#     
#     df = np.sqrt(cov_scale * df2)
#     y = model(x, *p)
#     delta = t_val * df
#     upperband = y + delta
#     lowerband = y- delta
#     
#     return [y, upperband, lowerband]
#             
# dfdp = [lambda a, x, b :x**b, lambda a, x , b: a * x**b * np.log(x)] 
# 
# nom = unp.nominal_values(py)
# std = unp.std_devs(py)
# 
# px = np.linspace(2,6,5)
# # plot the nominal value
# plt.plot(px, nom, c='r')
# 
# print(np.mean(np.square(nom - 2 * std)))
# print(np.mean(np.square(nom + 2 * std)))
# # And the 2sigma uncertaintie lines
# plt.plot(px, nom - 2 * std, c='c')
# plt.plot(px, nom + 2 * std, c='c')
# band, upper, lower = confidence_band(x, dfdp, 0.05, func)

