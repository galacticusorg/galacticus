'use strict';

window.Validation = (function () {

  function renderHeader(data) {
    document.getElementById('last-update').textContent = new Date(data.commit.timestamp).toString();

    const repoLink = document.getElementById('repository-link');
    repoLink.href = data.repoUrl;
    repoLink.textContent = data.repoUrl;

    const commitLink = document.getElementById('commit-link');
    commitLink.href = data.commit.url;
    commitLink.textContent = data.commit.id;

    const parameterLink = document.getElementById('parameter-link');
    if (parameterLink) {
      parameterLink.href = data.repoUrl + '/blob/' + data.commit.id + '/' + data.parameterFile;
      parameterLink.textContent = data.parameterFile;
    }

    const commandLine = document.getElementById('command-line');
    if (commandLine) {
      commandLine.textContent = data.commandLine;
    }
  }

  function attachDownloadButton(data, filename) {
    const btn = document.getElementById('dl-button');
    if (!btn) return;
    btn.onclick = () => {
      const a = document.createElement('a');
      a.href = 'data:,' + JSON.stringify(data, null, 2);
      a.download = filename || 'analyses_data.json';
      a.click();
    };
  }

  function renderGraph(parent, dataset) {
    const canvas = document.createElement('canvas');
    canvas.className = 'benchmark-chart';
    parent.appendChild(canvas);

    const src = dataset.data;
    const yIsLog = dataset.attributes.yAxisIsLog == 1;
    const xIsLog = dataset.attributes.xAxisIsLog == 1;

    const data = [];
    const dataTarget = [];
    let yMin =  Infinity;
    let yMax = -Infinity;

    for (let i = 0; i < src.xDataset.length; i++) {
      if (src.yDataset[i] != 0.0) {
        const p = { x: src.xDataset[i], y: src.yDataset[i] };
        if ('yError' in src) {
          p.yMin = src.yDataset[i] - src.yError[i];
          p.yMax = src.yDataset[i] + src.yError[i];
        }
        data.push(p);
        if (p.yMin < yMin && (!yIsLog || p.yMin > 0.0)) yMin = p.yMin;
        if (p.yMax > yMax) yMax = p.yMax;
      }
      if (src.yDatasetTarget[i] != 0.0) {
        const p = { x: src.xDataset[i], y: src.yDatasetTarget[i] };
        if ('yErrorTarget' in src) {
          p.yMin = src.yDatasetTarget[i] - src.yErrorTarget[i];
          p.yMax = src.yDatasetTarget[i] + src.yErrorTarget[i];
        }
        dataTarget.push(p);
        if (p.yMin < yMin && (!yIsLog || p.yMin > 0.0)) yMin = p.yMin;
        if (p.yMax > yMax) yMax = p.yMax;
      }
    }

    const options = {
      scales: {
        x: {
          ticks: {
            beginAtZero: false,
            callback: (val) => val.toExponential(2),
            maxRotation: 45,
            minRotation: 45
          },
          title: { display: true, text: dataset.attributes.xAxisLabel }
        },
        y: {
          ticks: {
            beginAtZero: false,
            callback: (val) => val.toExponential(2)
          },
          title: { display: true, text: dataset.attributes.yAxisLabel },
          min: yMin,
          max: yMax
        }
      }
    };
    if (xIsLog) options.scales.x.type = 'logarithmic';
    if (yIsLog) options.scales.y.type = 'logarithmic';

    new Chart(canvas, {
      type: 'scatterWithErrorBars',
      data: {
        datasets: [
          {
            label: 'Galacticus',
            data: data,
            borderColor: '#ff0000',
            backgroundColor: '#ffff00'
          },
          {
            label: dataset.attributes.targetLabel,
            data: dataTarget,
            borderColor: '#483D8B',
            backgroundColor: '#6495ED'
          }
        ]
      },
      options: options
    });
  }

  function renderAnalyses(dataSets, parent) {
    const main = parent || document.getElementById('main');
    for (const dataSet of dataSets) {
      const setElem = document.createElement('div');
      setElem.className = 'analysis';
      main.appendChild(setElem);

      const nameElem = document.createElement('h1');
      nameElem.className = 'analysis-title';
      nameElem.textContent = dataSet.attributes.name;
      setElem.appendChild(nameElem);

      const graphsElem = document.createElement('div');
      graphsElem.className = 'analysis-graphs';
      setElem.appendChild(graphsElem);

      renderGraph(graphsElem, dataSet);
    }
  }

  function renderStandard(data) {
    renderHeader(data);
    attachDownloadButton(data);
    renderAnalyses(data.results);
  }

  return {
    renderHeader,
    attachDownloadButton,
    renderGraph,
    renderAnalyses,
    renderStandard
  };
})();
