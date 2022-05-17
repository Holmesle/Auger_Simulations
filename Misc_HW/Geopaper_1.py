
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

labels = 'Confirmed Physical Hazard', 'Unconfirmed Physical Hazard', 'Unconfirmed Environmental Hazard', 'Confirmed Environmental Hazard', 'Mitigated or Non-hazardous'
sizes = [6439, 60279, 21262, 1363, 51309]
print(sum(sizes))
c = ['coral', 'lightseagreen', 'darkturquoise', 'salmon', 'dodgerblue']
explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

# fig1 = plt.figure(figsize=[5,5])
# plt.rcParams["font.family"] = "Times New Roman"
# plt.pie(sizes, colors = c,labels = labels, autopct='%1.1f%%', startangle=90)
# plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

file = r'C:\Users\holme\OneDrive\Desktop\python_code\Misc_HW\table.csv'

with open(file, 'r') as myfile:
    d = dict(pd.read_csv(myfile, dtype = float, delimiter=','))
    Year = d['Sector']
    Coal = d['Coal']
    Metal = d['Metal']
    Nonmetal = d['Nonmetal']
    Stone = d['Stone']
    Sand_n_Gravel = d['Sand & gravel']


def first_derivative(x, f):  # forward difference approximation
    h = x[0] - x[1]
    f_prime = np.zeros(len(f))
    for i in range(1,len(f)):
        numerator = f[i] - f[i-1] 
        f_prime[i] = numerator/(h)
    return f_prime


Total = Coal + Metal + Nonmetal + Sand_n_Gravel
Type = [Coal, Metal, Nonmetal, Sand_n_Gravel, Total]
Label = ['Coal', 'Metal','Nonmetal','Sand and Gravel', 'Total']
# ROC_coal = central_dif(Coal, np.arange(0, len(Year)))
# plt.plot(Year, central_dif(Coal,ROC_coal))
t = np.arange(0,len(Coal),1)
i = -1
ROC_T = first_derivative(t, Total)
plt.plot(Year, ROC_T, 'o', color=c[i], label=Label[i])
plt.plot(Year, ROC_T, '--', color=c[i], label = Label[i])
plt.xlabel('Year')
plt.ylabel('ROC in Number of Active Mines')
# plt.figure()
# for i in range(len(Type)):
#     plt.plot(Year,Type[i],'--',color = c[i])
#     plt.plot(Year, Type[i], 'o', color=c[i], label = Label[i])
# plt.legend()
# plt.xlabel('Year')
# plt.ylabel('Number of Active Mines')
# plt.show()
print(Total)
print('l')
