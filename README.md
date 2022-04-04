# SortStream: Sort-aware Gaussian-weighted Streamgraph
Here is the source code for the SortStream

You can find the layout algorithm in Layout_Ours.js

The input json file should be formatted as below:
```
[
    {
        "name": "a",
        "fillcolor": "rgb(10,20,30)",
        "size": [1, 2 ,3 ,4 ,5 ,6 ] // All data should be of the same length
    },
    {
        "name": "b",
        "fillcolor": "rgb(40,50,60)", 
        "size": [0,0,1,2,3,0] // Data of different lengths needs to be filled with zeros
    }
]
```
This paper was submitted to IEEE Transactions on Visualization and Computer Graphics 2022
### Abstract
Streamgraphs have been widely used in information visualization because of their beautiful display of multiple time
series. However, how to improve the readability of streamgraph is a long-standing research problem due to the influence of
data stacking. In this study, we propose a new technique to create streamgraphs, which is called SortStream. To make the
streamgraphs smoother, we extend the traditional method of assigning weights based on data thickness to a sort-aware Gaussian
weight assignment method to make the middle layer flatter and improve the readability of the overall graph. In addition, in terms
of layer sorting, we divide the position of layer fluctuation, carry out weighted calculations for baseline compensation between
layers, obtain layer sorting through a hierarchical clustering algorithm, and dynamically mine the weight of the relative length
and thickness of all layers in the clustering process to obtain better sorting. We provide two quantitative evaluations and user
experiments demonstrating that SortStream improves the readability and aesthetics of stream diagrams.





