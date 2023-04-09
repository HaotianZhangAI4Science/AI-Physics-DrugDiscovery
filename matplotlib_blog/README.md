[(60条消息) [seaborn\] seaborn学习笔记3-直方图Histogramplot_seaborn histogram_落痕的寒假的博客-CSDN博客](https://blog.csdn.net/LuohenYJ/article/details/90704424)

#### stick related 

```python
# set the stick interval as 2
from matplotlib.pyplot import MultipleLocator
x_major_locator=MultipleLocator(2)
ax.xaxis.set_major_locator(x_major_locator)

# set the stick size 
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# set the sticj range
ax.set_xlim([0,14])
```



