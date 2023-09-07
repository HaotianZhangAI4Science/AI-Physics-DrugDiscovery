[(60条消息) [seaborn\] seaborn学习笔记3-直方图Histogramplot_seaborn histogram_落痕的寒假的博客-CSDN博客](https://blog.csdn.net/LuohenYJ/article/details/90704424)



#### hist plot

```python
# the classical histgram plot according to the given intervals 
bins = np.arange(0,110,10)
plt.hist(all_per, bins, edgecolor='black')
plt.show()
```



#### stick related 

```python
# set the stick interval as 2
from matplotlib.pyplot import MultipleLocator
x_major_locator=MultipleLocator(2)
ax.xaxis.set_major_locator(x_major_locator)

# set the stick size 
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# set the stick range
ax.set_xlim([0,14])

# set the stick size and width
tick_size = 10
ax.tick_params(size=tick_size, width=3)
```

#### spines related 

```python
# set the color of spines as black 
for spine in ax.spines.values():
    spine.set_color('black')

# set spines invisible
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# or 
for spine in ax.spines.values():
    spine.set_visible(False)
# Even if you remove the spines, but the ticks, the number, still pertain. 
ax.set_xticks([])
ax.set_yticks([])
# set line width of spines
width = 3
for spine in ax.spines.values():
    spine.set_linewidth(width)
```

#### plot related

```python
# bar plot
counts = np.bincount(prop)/len(prop)
ax.bar(range(13), counts, width=0.6, align='center', color='white', ec=sns.color_palette()[0], linewidth=2)
# width is the width of bar
# ec is the color of outline of bar
# linewidth is the the line width of bar
```

#### scatter plot

```python
# using the label to directly assign the legend belonging. 
ax.scatter(actual_values, predicted_values, label='Predicted vs Actual')
```

#### Contour Plot

```python
sns.kdeplot(moses_transformed_new[:,0], moses_transformed_new[:,1], cmap="Blues", shade=True, shade_lowest=False, ax=ax, label='Moses')
```



#### other 

```python
# delete the grid in the ax
ax.grid(False)
```



##### Cmap (coloar-map) options 

```python
'viridis': Perceptually uniform colormap, suitable for accurately representing data in heatmaps.
'plasma': Another perceptually uniform colormap, but with a more vibrant color range.
'inferno': Similar to 'plasma', but with darker low values.
'magma': A softer, less contrasty perceptually uniform colormap.
'cividis': Designed to be legible by those with color vision deficiencies.
'Greys': A grayscale colormap.
'Blues': Varying shades of blue.
'Reds': Varying shades of red.
'coolwarm': Blue to white to red transition, good for diverging data.
'seismic': Blue to white to red, another option for diverging data.
'rainbow': Covers the spectrum of colors, but not perceptually uniform.
'jet': A popular colormap but often criticized for not being perceptually uniform.
```

