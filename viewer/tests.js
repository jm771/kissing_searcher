const testResults = document.getElementById('test-results');

function assert(condition, message) {
    const li = document.createElement('li');
    if (condition) {
        li.className = 'pass';
        li.textContent = `PASS: ${message}`;
    } else {
        li.className = 'fail';
        li.textContent = `FAIL: ${message}`;
    }
    testResults.appendChild(li);
}

function assertArraysAlmostEqual(a, b, epsilon = 1e-9) {
    if (a.length !== b.length) {
        return false;
    }
    for (let i = 0; i < a.length; i++) {
        if (Array.isArray(a[i]) && Array.isArray(b[i])) {
            if (!assertArraysAlmostEqual(a[i], b[i], epsilon)) {
                return false;
            }
        } else if (Math.abs(a[i] - b[i]) > epsilon) {
            return false;
        }
    }
    return true;
}

// --- Tests ---

assert(assertArraysAlmostEqual(identityMatrix(3), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 'identityMatrix(3)');

const m1 = [[1, 2], [3, 4]];
const m2 = [[5, 6], [7, 8]];
assert(assertArraysAlmostEqual(matrixMultiply(m1, m2), [[19, 22], [43, 50]]), 'matrixMultiply');

const v = [1, 2, 3];
const m = [[1, 0, 0], [0, 2, 0], [0, 0, 3]];
assert(assertArraysAlmostEqual(applyMatrix(m, v), [1, 4, 9]), 'applyMatrix');

const rot = rotationMatrix(3, 0, 1, Math.PI / 2);
assert(assertArraysAlmostEqual(rot, [[0, -1, 0], [1, 0, 0], [0, 0, 1]]), 'rotationMatrix');

const a = [1, 2, 3];
const b = [4, 5, 6];
assert(dot(a, b) === 32, 'dot');
assert(assertArraysAlmostEqual(add(a, b), [5, 7, 9]), 'add');
assert(assertArraysAlmostEqual(sub(a, b), [-3, -3, -3]), 'sub');
assert(assertArraysAlmostEqual(scale(a, 2), [2, 4, 6]), 'scale');

const p = [1, 2, 3];
const v1 = [1, 0, 0];
const w1 = [0, 1, 0];
const cam = [0, 0, 2];
const [x, y, s] = projectTo2D(p, v1, w1, cam);
assert(Math.abs(x - 1) < 1e-9 && Math.abs(y - 2) < 1e-9, 'projectTo2D');
