import torch

torch.argmin()
返回扁平化张量的最小值或沿某一维度的最小值的索引
这是 torch.min() 返回的第二个值。有关此方法的确切语义，请参见其文档。

torch.multinomial()
进行权重采样的函数
weights = torch.Tensor([0.2, 0.2, 0.3, 0.3])   # 采样权重
torch.multinomial(weights, 2)
