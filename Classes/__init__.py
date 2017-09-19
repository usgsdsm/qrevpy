# class Adam(object):
#     
#     def __init__(self):
#         self.test2 = None
# 
# 
# class Ben(object):
#     
#     def __init__(self):
#         self.test = Adam()
#         
#         
# a = Ben()
# 
# ref = getattr(a, 'test')
# ref2 = getattr(ref, 'test2')
# 
# ref2 = 'a'
# 
# print(a.test.test2)

# import numpy as np
# 
# a = np.array([2,2.5,3,2.7])
# b = [2,3,4,6]
# c = [1,2,3,4,5,6,7, 8]
# 
# d = np.interp(c,b,a)
# print (d)