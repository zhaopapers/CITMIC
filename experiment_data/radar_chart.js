//TCGA
option = {
  title: {
    text: 'TCGA'
  },
  tooltip: {
    trigger: 'axis'
  },
  legend: {
    data: ['CT-TME score', 'CT-TME score III/IV','ssGSEA','ssGSEA III/IV']
  },
  grid: {
    left: '3%',
    right: '4%',
    bottom: '3%',
    containLabel: true
  },
  toolbox: {
    feature: {
      saveAsImage: {}
    }
  },
  xAxis: {
    type: 'category',
    boundaryGap: false,
    data: ['1', '2', '3', '4', '5', '6', '7','8','9','10']
  },
  yAxis: {
    type: 'value'
  },
  series: [
 
    {
      name: 'CT-TME score',
      type: 'line',
    
        
      data: [0.713,0.708, 0.690, 0.689, 0.713 , 0.685 , 0.703,0.706 ,0.662 ,0.665],
      type: 'line',
      symbol: 'circle',
      symbolSize: 15,   
           lineStyle: {
        color: '#001852',
        width: 2,
        type: 'solid'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#001852',
        color: 'white'
      }
    },
 
    {
      name: 'CT-TME score III/IV',
      type: 'line',

      data: [0.767,0.761,0.782,0.810,0.811,0.809,0.845,0.865,0.794,0.830],
      
      symbol: 'circle',
      symbolSize: 15,   
           lineStyle: {
        color: '#e01f54',
        width:2,
        type: 'solid'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#e01f54',
        color: 'white'
      }
    },
    
    
    
        {
      name: 'ssGSEA',
      type: 'line',
  
        
      data: [0.731,0.665, 0.625, 0.635, 0.641 , 0.623 , 0.632,0.634 ,0.610 ,0.636],
      type: 'line',
      symbol: 'circle',
      symbolSize: 15,   
           lineStyle: {
        color: '#001852',
        width: 2,
        type: 'dashed'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#001852',
        color: 'white'
      }
    },
  
    {
      name: 'ssGSEA III/IV',
      type: 'line',
  
      data: [0.720,0.670,0.652,0.639,0.622,0.641,0.654,0.665,0.671,0.700],
      
      symbol: 'circle',
      symbolSize: 15,   
           lineStyle: {
        color: '#e01f54',
        width: 2,
        type: 'dashed'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#e01f54',
        color: 'white'
      }
    },
    
  ]
};

//GSE22155
  option = {
  title: {
    text: 'GSE22155'
  },
  tooltip: {
    trigger: 'axis'
  },
  legend: {
    data: ['CT-TME score III/IV','ssGSEA III/IV']
  },
  grid: {
    left: '3%',
    right: '4%',
    bottom: '3%',
    containLabel: true
  },
  toolbox: {
    feature: {
      saveAsImage: {}
    }
  },
  xAxis: {
    type: 'category',
    boundaryGap: false,
    data: ['1', '2', '3', '4']
  },
  yAxis: {
    type: 'value'
  },
  series: [
 
 
    {
      name: 'CT-TME score III/IV',
      type: 'line',

      data: [0.755,0.716,0.567, 0.631],
      
      symbol: 'circle',
      symbolSize: 15,   
           lineStyle: {
        color: '#e01f54',
        width:2,
        type: 'solid'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#e01f54',
        color: 'white'
      }
    },
    
    
    
       
    {
      name: 'ssGSEA III/IV',
      type: 'line',
  
      data: [0.507,0.604,0.508,0.566],
      
      symbol: 'circle',
      symbolSize: 15,   
      lineStyle: {
        color: '#e01f54',
        width: 2,
        type: 'dashed'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#e01f54',
        color: 'white'
      }
    },
    
  ]
};

//GSE19234
  option = {
  title: {
    text: 'GSE19234'
  },
  tooltip: {
    trigger: 'axis'
  },
  legend: {
    data: ['CT-TME score III/IV','ssGSEA III/IV']
  },
  grid: {
    left: '3%',
    right: '4%',
    bottom: '3%',
    containLabel: true
  },
  toolbox: {
    feature: {
      saveAsImage: {}
    }
  },
  xAxis: {
    type: 'category',
    boundaryGap: false,
    data: ['1', '2', '3', '4']
  },
  yAxis: {
    type: 'value'
  },
  series: [
 
 
    {
      name: 'CT-TME score III/IV',
      type: 'line',

      data: [0.8333333,0.5820431,0.6282323,0.6425606,0.7720148,0.7720148,0.7928478,0.8299473,0.8299473],
      
      symbol: 'circle',
      symbolSize: 15,   
           lineStyle: {
        color: '#e01f54',
        width:2,
        type: 'solid'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#e01f54',
        color: 'white'
      }
    },
    
    
    
       
    {
      name: 'ssGSEA III/IV',
      type: 'line',
  
      data: [0.8095238,0.4157855,0.5782555,0.5321772,0.4892586,0.4892586,0.5039504,0.4681315,0.4681315],
      
      symbol: 'circle',
      symbolSize: 15,   
      lineStyle: {
        color: '#e01f54',
        width: 2,
        type: 'dashed'
      },
      itemStyle: {
        borderWidth: 3,
        borderColor: '#e01f54',
        color: 'white'
      }
    },
    
  ]
};
