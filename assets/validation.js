'use strict';

window.Validation = (function () {

  const charts = [];

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

  function setupControls() {
    const main = document.getElementById('main');
    if (!main) return;

    const controls = document.createElement('div');
    controls.className = 'plot-controls';
    controls.style.cssText = 'padding: 15px; background: #f4f4f4; border-radius: 8px; margin-bottom: 20px; display: flex; gap: 30px; align-items: center; border: 1px solid #ddd;';
    
    controls.innerHTML = `
      <div style="display: flex; flex-direction: column; gap: 5px;">
        <label style="font-weight: bold; font-size: 0.9em;">Point Size</label>
        <input type="range" id="point-size-slider" min="1" max="15" value="4" style="cursor: pointer;">
      </div>
      <div style="display: flex; flex-direction: column; gap: 5px;">
        <label style="font-weight: bold; font-size: 0.9em;">Line/Error Thickness</label>
        <input type="range" id="line-thickness-slider" min="1" max="8" value="2" style="cursor: pointer;">
      </div>
    `;
    main.prepend(controls);

    document.getElementById('point-size-slider').addEventListener('input', (e) => {
      const val = parseInt(e.target.value);
      charts.forEach(chart => {
        chart.data.datasets.forEach(ds => {
          ds.pointRadius = val;
          ds.radius = val;
          ds.pointHoverRadius = val + 2;
        });
        chart.update('none');
      });
    });

    document.getElementById('line-thickness-slider').addEventListener('input', (e) => {
      const val = parseInt(e.target.value);
      charts.forEach(chart => {
        chart.data.datasets.forEach(ds => ds.borderWidth = val);
        chart.update('none');
      });
    });
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
      responsive: true,
      maintainAspectRatio: true,
      aspectRatio: 1.5,
      plugins: {
        legend: {
          position: 'top',
          labels: {
            usePointStyle: true,
            padding: 20,
            font: { size: 12 }
          }
        }
      },
      scales: {
        x: {
          grid: { color: 'rgba(0, 0, 0, 0.05)' },
          ticks: {
            beginAtZero: false,
            callback: (val) => {
              const absVal = Math.abs(val);
              if (absVal === 0) return '0';
              if (absVal < 1e-3 || absVal >= 1e4) return val.toExponential(1);
              return val.toLocaleString();
            }
          },
          title: { display: true, text: dataset.attributes.xAxisLabel, font: { weight: 'bold' } }
        },
        y: {
          grid: { color: 'rgba(0, 0, 0, 0.05)' },
          ticks: {
            beginAtZero: false,
            callback: (val) => {
              const absVal = Math.abs(val);
              if (absVal === 0) return '0';
              if (absVal < 1e-3 || absVal >= 1e4) return val.toExponential(1);
              return val.toLocaleString();
            }
          },
          title: { display: true, text: dataset.attributes.yAxisLabel, font: { weight: 'bold' } },
          min: yMin,
          max: yMax
        }
      }
    };
    if (xIsLog) options.scales.x.type = 'logarithmic';
    if (yIsLog) options.scales.y.type = 'logarithmic';

    const chart = new Chart(canvas, {
      type: 'scatterWithErrorBars',
      data: {
        datasets: [
          {
            label: 'Galacticus',
            data: data,
            borderColor: '#0072B2',
            backgroundColor: 'rgba(0, 114, 178, 0.3)',
            pointStyle: 'circle',
            radius: 4,
            pointRadius: 4,
            borderWidth: 2
          },
          {
            label: dataset.attributes.targetLabel,
            data: dataTarget,
            borderColor: '#000000',
            backgroundColor: 'rgba(0, 0, 0, 0.1)',
            pointStyle: 'rectRot',
            radius: 5,
            pointRadius: 5,
            borderWidth: 2
          }
        ]
      },
      options: options
    });
    charts.push(chart);
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
    setupControls();
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
