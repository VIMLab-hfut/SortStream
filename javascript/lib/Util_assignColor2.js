// function assignColor(layers) {
//     getAdjacentMatrix_StreamGraph(layers)
//     // layers.forEach((d,i)=>d.fillcolor = "red")
//     return layers
// }



function assignLayersColor_4Color(layers, layersData) {
    let Palettes = ['rgb(247,252,253)', 'rgb(229,245,249)', 'rgb(204,236,230)', 'rgb(153,216,201)', 'rgb(102,194,164)', 'rgb(65,174,118)', 'rgb(35,139,69)', 'rgb(0,109,44)', 'rgb(0,68,27)']
    // let Palettes = ['rgb(166,206,227)', 'rgb(31,120,180)', 'rgb(178,223,138)', 'rgb(51,160,44)', 'rgb(251,154,153)', 'rgb(227,26,28)', 'rgb(253,191,111)', 'rgb(255,127,0)', 'rgb(202,178,214)', 'rgb(106,61,154)', 'rgb(255,255,153)', 'rgb(177,89,40)']
    //     'rgb(105,105,133)', 'rgb(154,156,204)', 'rgb(138,138,174)', 'rgb(152,151,195)',
    //     'rgb(120,117,146)', 'rgb(150,149,181)', 'rgb(129,127,164)', 'rgb(127,124,153)',
    //     'rgb(163,164,208)', 'rgb(137,140,173)', 'rgb(122,121,153)', 'rgb(136,135,169)',
    //     'rgb(105,102,123)', 'rgb(103,100,122)',  'rgb(87,86,104)',
    //     'rgb(91,92,112)', 'rgb(157,159,208)', 'rgb(85,88,97)', 'rgb(109,107,129)',
    //     'rgb(126,126,160)', 'rgb(150,149,191)', 'rgb(92,94,117)', 'rgb(169,168,212)',
    //     'rgb(173,172,216)', 'rgb(97,95,117)', 'rgb(162,163,207)', 'rgb(135,137,175)',
    //     'rgb(132,132,170)', 'rgb(110,111,131)'
    // ]

    Palettes = shuffle(Palettes)
    let thisPalette = Palettes

    let adjacencyMatrix = getAdjacentMatrix_StreamGraph(layers)

    // let count_colors = getMax_Color_Num(adjacencyMatrix)
    let count_colors = Palettes.length

    let this_color = getRandom(1, count_colors),
        this_color_index = 1,
        this_layer_index = 1,
        past_layer_index = 0;

    let array_colors = getArray2D(layers.length, 2, "a");
    array_colors[0][0] = 1;
    array_colors[0][1] = 1;

    while (this_layer_index < array_colors.length) {
        while (this_color_index <= count_colors) {
            if (this_layer_index >= array_colors.length) {
                break
            }
            past_layer_index = 0
            while (past_layer_index < this_layer_index) {
                if (adjacencyMatrix[this_layer_index][past_layer_index] === 0) {
                    past_layer_index++
                    continue
                } else {
                    if (array_colors[past_layer_index][0] === this_color) {
                        break
                    } else {
                        past_layer_index++
                        continue
                    }
                }
            }
            if (past_layer_index < this_layer_index) {
                // this_color++
                this_color = this_color + 1 > count_colors ? 1 : this_color + 1
                this_color_index++
            } else {
                array_colors[this_layer_index][0] = this_color
                array_colors[this_layer_index][1] = this_color_index
                this_layer_index++
                // this_color = 1
                this_color = this_color + 1 > count_colors ? 1 : this_color + 1
                this_color_index = 1
            }
        }
        if (this_color_index > count_colors) {
            this_layer_index--
            this_color = array_colors[this_layer_index][0] + 1 > count_colors ? 1 : array_colors[this_layer_index][0] + 1
            this_color_index = array_colors[this_layer_index][1] + 1
        }
    }

    for (let i = 0; i < array_colors.length; i++) {
        layersData[i].fillcolor = thisPalette[array_colors[i][0] - 1]
        layers[i].fillcolor = thisPalette[array_colors[i][0] - 1]
    }

    return layersData
}


/**
 * ????????????streamgraph???????????????
 * @param {?????????} layers
 */
function getAdjacentMatrix_StreamGraph(layers) {
    let resMatrix = getArray2D(layers.length, layers.length, 0)
    for (let i = 0; i < layers.length - 1; i++) { //????????????layer
        for (let j = 1; j < layers[i].size.length; j++) {
            if (layers[i].size[j - 1] + layers[i].size[j] === 0) {
                continue
            } else {
                for (let k = i + 1; k < layers.length; k++) {
                    if (layers[k].size[j - 1] + layers[k].size[j] > 0) {
                        // resMatrix[i][k] += getPoint_Distance_Graph(j - 1, layers[k].yBottom[j - 1], j, layers[k].yBottom[j])
                        resMatrix[i][k] = 1
                        resMatrix[k][i] = resMatrix[i][k]
                        break
                    }
                }
            }
        }
    }
    return resMatrix;
}


function assignLayersColor_4Color2(layers) {


    let Palettes = ['rgb(166,206,227)', 'rgb(31,120,180)', 'rgb(178,223,138)', 'rgb(51,160,44)', 'rgb(251,154,153)', 'rgb(227,26,28)', 'rgb(253,191,111)', 'rgb(255,127,0)', 'rgb(202,178,214)', 'rgb(106,61,154)', 'rgb(255,255,153)', 'rgb(177,89,40)']
    // let Palettes = ['rgb(237,248,251)', 'rgb(191,211,230)', 'rgb(158,188,218)', 'rgb(140,150,198)', 'rgb(136,86,167)', 'rgb(129,15,124)']
    // let Palettes = ['rgb(247,247,247)', 'rgb(217,217,217)', 'rgb(189,189,189)', 'rgb(150,150,150)', 'rgb(99,99,99)', 'rgb(37,37,37)']
    // let Palettes = ['rgb(255,255,204)', 'rgb(217,240,163)', 'rgb(173,221,142)', 'rgb(120,198,121)', 'rgb(49,163,84)', 'rgb(0,104,55)']
    // let thisPalette = shuffle(Palettes).slice(0, 4)
    // let thisPalette = shuffle(Palettes)
    // let thisPalette =
    // let thisPalette = ['rgb(247,247,247)', 'rgb(204,204,204)', 'rgb(150,150,150)', 'rgb(99,99,99)', 'rgb(37,37,37)'].slice(0, 4)
    // let thisPalette = Palettes.slice(Math.round(Palettes.length / 4), Math.round(Palettes.length / 4) + 4)

    for (let i = 0; i < layers.length; i++) {
        let color = getRandom(0, Palettes.length - 1);
        layers[i].fillcolor = Palettes[color]
        // Palettes.slice(color,1)
    }
    return layers
}

function getMax_Color_Num(adjacencyMatrix) {
    let res = 4
    for (let i = 0; i < adjacencyMatrix.length; i++) {
        let connectedPoints = A_Copy_Of(['a'])
        connectedPoints.shift()
        for (let j = 0; j < i; j++) {
            if (adjacencyMatrix[i][j] === 1 && i !== j) {
                connectedPoints.push(j)
            }
        }
        if (connectedPoints.length >= res) { //?????????4???????????????
            for (let k = res; k < connectedPoints.length + 1; k++) {
                let all_child_array = choose(connectedPoints, k)
                for (let l = 0; l < all_child_array.length; l++) {
                    all_child_array[l].push(i)
                    if (isCompleteSubGraph(all_child_array[l], adjacencyMatrix)) {
                        res = all_child_array[l].length;
                        break
                    }
                }
            }

        }
    }
    return res
}

function isCompleteSubGraph(points, adjacencyMatrix) {
    for (let i = 0; i < points.length - 1; i++) {
        for (let j = i + 1; j < points.length; j++) {
            if (adjacencyMatrix[points[i]][points[j]] === 1) {
                continue
            } else {
                return false
            }
        }
    }
    console.log(points);
    return true
}

/**
 * ??????????????????
 * @param {??????} arr
 * @param {????????????????????????} size
 */
function choose(arr, size) {
    var allResult = [];

    (function (arr, size, result) {
        var arrLen = arr.length;
        if (size > arrLen) {
            return;
        }
        if (size == arrLen) {
            allResult.push([].concat(result, arr))
        } else {
            for (var i = 0; i < arrLen; i++) {
                var newResult = [].concat(result);
                newResult.push(arr[i]);

                if (size == 1) {
                    allResult.push(newResult);

                } else {
                    var newArr = [].concat(arr);
                    newArr.splice(0, i + 1);
                    arguments.callee(newArr, size - 1, newResult);

                }

            }
        }
    })(arr, size, []);
    return allResult;
}