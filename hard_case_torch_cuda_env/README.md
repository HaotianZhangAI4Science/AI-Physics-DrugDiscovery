### ImportError: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.26‘ not found

```shell
# find whether this so.6 in your conda env
strings /home/app/anaconda3/lib/libstdc++.so.6 | grep GLIBCXX
# if so
export LD_LIBRARY_PATH=/home/app/anaconda3/lib
```

### CUDA error: CUBLAS STATUS_INVALID VALUE when calling cublasGemmEx( handle ...

```shell
# Once I tried to use hugging-face transformer and pytorch-lightning to train something, this bug occurred again and again. Meanwhile, the model was deseperated from the cuda device once I straded training. Then, I found the problem is that my torch version is not compatible with the cuda version. Annoying, right? 
```

