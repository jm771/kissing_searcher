// ===================================
// Core Mathematical Functions
// ===================================

function add(a, b) { return a.map((ai, i) => ai + b[i]); }
function sub(a, b) { return a.map((ai, i) => ai - b[i]); }
function scaleVector(a, s) { return a.map(ai => ai * s); }
function dot(a, b) { return a.reduce((sum, ai, i) => sum + ai * b[i], 0); }

function identityMatrix(n) {
  return Array(n).fill(null).map((_, i) =>
    Array(n).fill(null).map((_, j) => (i === j ? 1 : 0))
  );
}

function matrixMultiply(a, b) {
  const n = a.length;
  const result = Array(n).fill(null).map(() => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      for (let k = 0; k < n; k++) {
        result[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return result;
}

function applyMatrix(matrix, vector) {
  const n = vector.length;
  const result = Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result[i] += matrix[i][j] * vector[j];
    }
  }
  return result;
}

function rotationMatrix(n, axis1, axis2, theta) {
  const matrix = identityMatrix(n);
  const c = Math.cos(theta);
  const s = Math.sin(theta);
  matrix[axis1][axis1] = c;
  matrix[axis1][axis2] = -s;
  matrix[axis2][axis1] = s;
  matrix[axis2][axis2] = c;
  return matrix;
}

function projectTo2D(P, v, w, camera) {
  const x = dot(P, v);
  const y = dot(P, w);
  let resid = sub(P, scaleVector(v, x));
  resid = sub(resid, scaleVector(w, y));
  const distance_vec = sub(resid, camera);
  const scale_factor = Math.sqrt(dot(camera, camera) / dot(distance_vec, distance_vec));
  return [x, y, scale_factor];
}

// ===================================
// Browser Application Logic
// ===================================

let IMPORT_RADIUS = 1;

function runApp(frames) {
  const canvas = document.getElementById('canvas');
  const ctx = canvas.getContext('2d');

  let n = frames[0][0].length;
  let sphereCenters = frames[0].map(p => p.map(coord => coord / IMPORT_RADIUS));

  const frameSlider = document.getElementById('frameSlider');
  const frameSliderValue = document.getElementById('frameSliderValue');
  frameSlider.max = frames.length - 1;
  frameSlider.addEventListener('input', () => {
    const frameIndex = parseInt(frameSlider.value, 10);
    sphereCenters = frames[frameIndex].map(p => p.map(coord => coord / IMPORT_RADIUS));
    frameSliderValue.textContent = frameIndex;
    draw();
  });

  const jumpBackwardBtn = document.getElementById('jumpBackwardBtn');
  const jumpForwardBtn = document.getElementById('jumpForwardBtn');
  const frameJumpInput = document.getElementById('frameJumpInput');

  function jumpFrames(direction) {
    const jumpValue = parseInt(frameJumpInput.value, 10);
    let currentFrame = parseInt(frameSlider.value, 10);
    let newFrame = currentFrame + (direction * jumpValue);

    newFrame = Math.max(0, Math.min(newFrame, frameSlider.max));

    if (newFrame !== currentFrame) {
        frameSlider.value = newFrame;
        const event = new Event('input', { bubbles: true });
        frameSlider.dispatchEvent(event);
    }
  }

  jumpBackwardBtn.addEventListener('click', () => jumpFrames(-1));
  jumpForwardBtn.addEventListener('click', () => jumpFrames(1));

  let zoom = 0.5;

  let circleSize = 1;
  let rotationStep = 0.31415926;
  let tumble = false;
  let scaleByDistance = false;
  let rotation = identityMatrix(n);

  function resetRotation() {
    rotation = identityMatrix(n);
    draw();
  }

  document.getElementById('resetButton').addEventListener('click', resetRotation);

  document.getElementById('zoomSlider').addEventListener('input', (e) => {
    zoom = parseFloat(e.target.value);
    document.getElementById('zoomSliderValue').textContent = zoom.toFixed(2);
    draw();
  });

  document.getElementById('circleSizeSlider').addEventListener('input', (e) => {
    circleSize = parseFloat(e.target.value);
    document.getElementById('circleSizeSliderValue').textContent = circleSize.toFixed(2);
    draw();
  });

  document.getElementById('rotationStepSlider').addEventListener('input', (e) => {
    rotationStep = parseFloat(e.target.value);
    document.getElementById('rotationStepSliderValue').textContent = rotationStep.toPrecision(4);
  });

  document.getElementById('tumbleToggle').addEventListener('change', (e) => {
    tumble = e.target.checked;
    if (tumble) tumbleLoop();
  });

  document.getElementById('distanceToggle').addEventListener('change', (e) => {
    scaleByDistance = e.target.checked;
    draw();
  });

  function createRotationButtons(axis1, axis2, container) {
    const buttonGroup = document.createElement('div');
    buttonGroup.className = 'button-group';
    const label = document.createElement('span');
    label.className = 'rotation-label';
    label.textContent = `(${axis1 + 1},${axis2 + 1})`;
    const leftButton = document.createElement('button');
    leftButton.textContent = '<';
    leftButton.addEventListener('click', () => {
      const rot = rotationMatrix(n, axis1, axis2, -rotationStep);
      rotation = matrixMultiply(rotation, rot);
      draw();
    });
    const rightButton = document.createElement('button');
    rightButton.textContent = '>';
    rightButton.addEventListener('click', () => {
      const rot = rotationMatrix(n, axis1, axis2, rotationStep);
      rotation = matrixMultiply(rotation, rot);
      draw();
    });
    buttonGroup.appendChild(label);
    buttonGroup.appendChild(leftButton);
    buttonGroup.appendChild(rightButton);
    container.appendChild(buttonGroup);
  }

  function setupControls() {
    const xAxisControls = document.getElementById('x-axis-controls');
    const yAxisControls = document.getElementById('y-axis-controls');
    const planeRotationControl = document.getElementById('plane-rotation-control');
    xAxisControls.innerHTML = '';
    yAxisControls.innerHTML = '';
    planeRotationControl.innerHTML = '';
    for (let i = 2; i < n; i++) {
      createRotationButtons(0, i, xAxisControls);
      createRotationButtons(1, i, yAxisControls);
    }
    createRotationButtons(0, 1, planeRotationControl);
  }

  function draw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.save();
    ctx.translate(canvas.width / 2, canvas.height / 2);
    ctx.scale(zoom, zoom);

    const basis_v = applyMatrix(rotation, [1, 0, ...Array(n - 2).fill(0)]);
    const basis_w = applyMatrix(rotation, [0, 1, ...Array(n - 2).fill(0)]);
    const camera = applyMatrix(rotation, [0, 0, -2, ...Array(n - 3).fill(0)]);

    let projected = sphereCenters.map(p => {
      let main = projectTo2D(p, basis_v, basis_w, camera);
      let overlap = sphereCenters.map(x => {let d = sub(x, p); let m = dot(d, d); return m == 0 ? 10 : m;} ).reduce((x, y) => Math.min(x, y));
      main.push(overlap);
      return main;
    });

    projected.push([0, 0, 1])

     projected.sort((a, b) => a[2] - b[2]);

    projected.forEach(([x, y, scale_factor, mindist]) => {
      ctx.beginPath();
      const radius = circleSize * 200 * (scaleByDistance ? scale_factor / 2 : 1);
      ctx.arc(x * 400, y * 400, radius, 0, 2 * Math.PI);

      // if (x == 0 && y == 0) { 
      //   ctx.fillStyle = `rgba(0, 255, 0, ${0.5 * scale_factor})`
      // } else 
      if (mindist < 0.99) {
        let overlapScale = Math.sqrt(mindist); // * mindist;
        overlapScale = (overlapScale - 0.9) * 10
        let r = 255 * (1-overlapScale);
        let g = 255 * overlapScale;
        ctx.fillStyle = `rgba(${r}, ${g}, 0, ${0.5 * scale_factor})`;
      } else {
        ctx.fillStyle = `rgba(0, 255, 0, ${0.5 * scale_factor})`;
      }  
      
      ctx.fill();
    });
    

    ctx.restore();
  }

  const TUMBLE_SUSTAIN = 75;
  let tumbleFrame = TUMBLE_SUSTAIN;
  let tumbleMatrix = null;

  function tumbleLoop() {
    if (!tumble) return;

    if (tumbleFrame == TUMBLE_SUSTAIN)
    {
      const axis1 = Math.floor(Math.random() * n);
      let axis2 = Math.floor(Math.random() * (n - 1));
      if (axis2 >= axis1) axis2++;
      const direction = Math.floor(Math.random() * 2) * 2 - 1;

      tumbleMatrix = rotationMatrix(n, axis1, axis2, rotationStep * direction / 20);
      tumbleFrame = 0;
    }
    else
    {
      tumbleFrame++;
    }
    
    rotation = matrixMultiply(tumbleMatrix, rotation);
    draw();
    requestAnimationFrame(tumbleLoop);
  }

  setupControls();
  draw();
}

// ===================================
// Entry Point
// ===================================

if (typeof module !== 'undefined' && module.exports) {
  // Export functions for testing
  module.exports = { add, sub, scaleVector, dot, identityMatrix, matrixMultiply, applyMatrix, rotationMatrix, projectTo2D };
} else {
  // Run the application in the browser
  document.addEventListener('DOMContentLoaded', () => {
    fetch('frames.json')
      .then(response => response.json())
      .then(runApp)
      .catch(error => console.error('Error loading or running application:', error));
  });
}
