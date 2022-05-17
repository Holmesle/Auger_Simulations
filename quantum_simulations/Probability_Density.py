import numpy as np
import matplotlib.pyplot as plt

# Variables
L = 1
N_bins = 100
N_links = 2
N_samples = 10000

# Single Link Generation
thetas = 2*np.pi*np.random.rand(N_samples)
x_positions = L*np.cos(thetas)
y_positions = L*np.sin(thetas)

# averages
avg_x = sum(x_positions)/N_samples
avg_x_sqrd = sum(x_positions**2)/N_samples

# probability density
bins = np.zeros(N_bins)
bins_range = np.linspace(-N_links*L, N_links*L,N_bins)
bin_width = 2*L/N_bins

for i in range(0,N_samples):
    # print(x_positions[i])
    bin_number = int(((x_positions[i] + N_links*L)/(2*N_links*L))*N_bins)
    # print(bin_number)
    bins[bin_number] = bins[bin_number] + 1

norm = bins/(N_samples)
PD = norm/bin_width
sqrt_avg_x_sqrd = np.sqrt(avg_x_sqrd)
print(sqrt_avg_x_sqrd)

# plot
# plt.figure()
# plt.plot(x_positions,y_positions,'ro')
# plt.gca().set_aspect('equal')

# plt.figure(figsize = [5,4])
# out = np.histogram(x_positions)
# plt.plot(bins_range, PD ,color = 'black')
# plt.axvline(sqrt_avg_x_sqrd,color = 'red')


# Multi Link Generation
def Multi_Link(N_links,N_samples):
    thetas = 2*np.pi*np.random.rand(N_samples, N_links)
    x_positions = L*np.cos(thetas)
    y_positions = L*np.sin(thetas)
    xf_positions = x_positions.sum(axis=1)
    yf_positions = y_positions.sum(axis=1)
    sqrd_x_avg = sum(xf_positions**2)/N_samples
    return xf_positions, yf_positions,sqrd_x_avg

# xf_positions,yf_positions = Multi_Link(N_links, N_samples)

# print(xf_positions)

# # averages
# avg_x = sum(xf_positions)/N_samples
# avg_x_sqrd = sum(xf_positions**2)/N_samples

# probability density
# bins = np.zeros(N_bins)
# bins_range = np.linspace(-N_links*L, N_links*L, N_bins)
# bin_width = 2*L/N_bins


# for i in range(0, N_samples):
#     # print(x_positions[i])
#     bin_number = int(((xf_positions[i] + N_links*L)/(2*N_links*L))*N_bins)
#     bins[bin_number] = bins[bin_number] + 1

N_max = 1000
step_size = 1
avg_x_sqrd = []
N_array = []
for N in range(1,N_max-2,step_size):
    x = Multi_Link(N, 1000)[2]
    # print(x)
    avg_x_sqrd.append(x)
    N_array.append(N)

# avg_x_sqrd = np.array(Multi_Link(N, 1000)[2])

print('x = ',np.shape(avg_x_sqrd))
plt.plot(N_array,avg_x_sqrd)
plt.xlabel('N')
plt.ylabel('<x^2>')
plt.show()

# plt.figure(figsize = [5,4])
# out = np.histogram(xf_positions)
# plt.plot(bins_range, PD ,color = 'black')
# plt.axvline(sqrt_avg_x_sqrd,color = 'red')
# plt.show()
