from tkinter import font
from sklearn.neighbors import KernelDensity
import seaborn as sns
def make_data(N, f=0.3, rseed=1):
    rand = np.random.RandomState(rseed)
    x = rand.randn(N)
    x[int(f * N):] += 5
    return x

x = make_data(1000)


font1 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 14,
}
# instantiate and fit the KDE model
kde = KernelDensity(bandwidth=1.0, kernel='gaussian')
x_d = np.linspace(-4, 8, 1000)
kde.fit(x[:, None])

# score_samples returns the log of the probability density
logprob = kde.score_samples(x_d[:, None])

index_start = 100
index_end = 500
ax=plt.gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
plt.fill_between(x_d[:index_end], np.exp(logprob[:index_end]), alpha=0.7)
bins=50
#plt.hist(x, density=True, bins=50, alpha=0.5)
plt.plot(x_d, np.exp(logprob),color=sns.color_palette()[0], linewidth=2)

plt.ylim(0, 0.2)
plt.legend(['KDE'],loc='upper left',prop=font1)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)


# 做2D histogram 的等高线图以及plot新采样的点上去，控制输出的colorbar，以及坐标轴的间距，坐标轴保留的小数点
import matplotlib.ticker as mtick
meta = np.array([tica_concatenated_0[:50000,0],tica_concatenated_0[:50000,1]])
fig,ax=plt.subplots(1,1)

alpha=0.8
eps = 1e-4
bins = (40,40)
xlim = (-2,2.5)
ylim = (-2.6,2.0)
ranges = (xlim, ylim)
H, xedges, yedges = np.histogram2d(meta[0], meta[1] ,bins=bins, range=ranges)
H += eps
H = -np.log(H)
x = np.linspace(xedges[0], xedges[-1],bins[0])
y = np.linspace(yedges[0], yedges[-1],bins[1])
X, Y = np.meshgrid(x, y)
plt.contourf(X,Y,H.T, cmap='Blues_r', alpha=0.7)

p2 = plt.scatter(tica_concatenated_0[50001:50250,0],tica_concatenated_0[50001:50250,1],color=sns.color_palette()[3], alpha=0.8)
p3 = plt.scatter(tica_concatenated_4[50001:50578,0],tica_concatenated_4[50001:50578,1], color=sns.color_palette()[1],alpha=0.8)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
x_major_locator=MultipleLocator(1.0)
y_major_locator = MultipleLocator(1.0)
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.xlim(xlim)
plt.ylim(ylim)
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
plt.legend([p2,p3],['SDEGen','CREST'], loc='upper right',fontsize=14)
plt.savefig('1_SDEGen.png', dpi=400)
