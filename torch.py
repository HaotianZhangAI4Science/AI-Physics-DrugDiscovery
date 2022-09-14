import torch

torch.argmin()
返回扁平化张量的最小值或沿某一维度的最小值的索引
这是 torch.min() 返回的第二个值。有关此方法的确切语义，请参见其文档。


torch.multinomial()
进行权重采样的函数
weights = torch.Tensor([0.2, 0.2, 0.3, 0.3])   # 采样权重
torch.multinomial(weights, 2)
返回两个采样的样本


关于用batch 训练pkt和lig如何用mask方便的取值，分别做图神经处理，最后concat到一张图上的技巧
从GraphBP项目当中获得的灵感
我们有一个rec_mask, 组成是[True, True, .., False, False...]即蛋白相应的位置是True, lig相应的位置是False
# test
a = torch.zeros(10, dtype=torch.float32)
rec_mask = torch.tensor([0,0,1,1,1,0,1,1,1,1])
rec_mask = rec_mask==1
a[~rec_mask] = torch.tensor([27,27,27], dtype=torch.float32)
a[rec_mask] = torch.tensor([34], dtype=torch.float32).repeat(7)
# model 
z = batch['atom_type']
zh = torch.zeros(z.shape[0],self.hidden_channels, dtype=torch.float32)  #self.hidden_channel是pkt和lig经过各自变换以后的维度
rec_nask = batch['rec_mask']
z[rec_mask]=0
rec_z = z[rec_mask]
lig_z = z[~rec_mask][:,0].to(torch.long)
rec_zh = self.pkt_embedding(rec_z)
lig_zh = self.lig_embedding(lig_z)
zh[rec_mask] = rec_zh
zh[~rec_mask] = lig_zh 


有时间看一下torch.view的机制
以下的代码是判断两个矩阵是否完全相同的
torch.all(torch.eq(test, z.view(-1,4)))



#处理离散变量
atomic_numbers = torch.LongTensor([6,7,8,9,15,16,17])  # C N O F P S Cl
element = data.ligand_element.view(-1, 1) == self.atomic_numbers.view(1, -1)   # (N_atoms, N_elements)
