// Physical parameters
let Lx;
let Ly;
let k;
let cp;
let rho;
let alpha;

// Numerical parameters
let nx;
let ny;
let dx;
let dy;
let eps;
let maxIt;
let BC1;
let BC2;
let BC3;
let BC4;

// Boundary conditions
function boundaryCondition1(z, k, nx, ny, dx, dy, bc1, type) {
   for (let i = 0; i < nx; i++) {
      if (type == "Dirichlet") {
         z[i][0] = bc1;
      }
      if (type == "Neumann") {
         z[i][0] = z[i][1] + k*bc1*dy;
      }
   }
   return z;
}

function boundaryCondition2(z, k, nx, ny, dx, dy, bc2, type) {
   for (let j = 0; j < ny; j++) {
      if (type == "Dirichlet") {
         z[0][j] = bc2;
      }
      if (type == "Neumann") {
         z[0][j] = z[1][j] + k*bc2*dx;
      }
   }
   return z;
}

function boundaryCondition3(z, k, nx, ny, dx, dy, bc3, type) {
   for (let j = 0; j < ny; j++) {
      if (type == "Dirichlet") {
         z[nx-1][j] = bc3;
      }
      if (type == "Neumann") {
         z[nx-1][j] = z[nx-2][j] + k*bc3*dx;
      }
   }
   return z;
}

function boundaryCondition4(z, k, nx, ny, dx, dy, bc4, type) {
   for (let i = 0; i < nx; i++) {
      if (type == "Dirichlet") {
         z[i][ny-1] = bc4;
      }
      if (type == "Neumann") {
         z[i][ny-1] = z[i][ny-2] + k*bc4*dy;
      }
   }
   return z;
}

function edgesBC(z, nx, ny) {
   z[0][0] = (z[0][1]+z[1][0])/2;
   z[0][ny-1] = (z[0][ny-2]+z[1][ny-1])/2;
   z[nx-1][0] = (z[nx-1][1]+z[nx-2][0])/2;
   z[nx-1][ny-1] = (z[nx-2][ny-1]+z[nx-1][ny-2])/2;
   return z;
}

// Check for convergence
function checkConvergence(z, z0, nx, ny, eps) {
   let maxDiff;
   for (let i = 1; i < nx - 1; i++) {
      for (let j = 1; j < ny - 1; j++) {
         maxDiff = Math.abs(z[i][j]-z0[i][j]);
         if (maxDiff < eps) {
            return true;
         }
      }
   }
   return false;
}

// Solve differential equation
function solve(z, z0, nx, ny, dx, dy, bc1, bc2, bc3, bc4, bc1type, bc2type, bc3type, bc4type, eps, maxIt) {

   for (let it = 0; it < maxIt; it++) {
      z0 = JSON.parse(JSON.stringify(z)); // deep copy
      for (let i = 1; i < nx-1; i++) {
         for (let j = 1; j < ny-1; j++) {
            z[i][j] = (dy*dy*(z[i+1][j]+z[i-1][j])+dx*dx*(z[i][j+1]+z[i][j-1]))/(2*(dx*dx+dy*dy));
         }
      }
      z = boundaryCondition1(z, k, nx, ny, dx, dy, bc1, bc1type);
      z = boundaryCondition2(z, k, nx, ny, dx, dy, bc2, bc2type);
      z = boundaryCondition3(z, k, nx, ny, dx, dy, bc3, bc3type);
      z = boundaryCondition4(z, k, nx, ny, dx, dy, bc4, bc4type);
      z = edgesBC(z, nx, ny);
      if (checkConvergence(z, z0, nx, ny, eps) == true) {
         alert("Solution converged after " + it + " iterations.")
         return z;
      }
   }
   alert("Solution did not converged!");
   return z;
}

function plot(Z) {
   let data = [{
      z: Z,
      zsmooth: 'best',
      type: 'heatmap',
      showscale: true,
      connectgaps: true,
      xaxis: 'x4',
      yaxis: 'y4',
      colorscale: 'Jet',
    }];
    
    let layout = {
      paper_bgcolor: '#121212',
      font_color: 'rgb(255,255,255)',
      //paper_bgcolor='rgb(233,233,233)'
      title: 'Contour plot',
      font: {
         family: 'sans-serif',
         size: 18,
         color: '#cccccc'
       }
    };
    
    Plotly.newPlot('plot', data, layout);
}

function compute() {
   Lx = +document.getElementById("Lx").value;
   Ly = +document.getElementById("Ly").value;
   k = +document.getElementById("k").value;
   cp = +document.getElementById("cp").value;
   rho = +document.getElementById("rho").value;
   alpha = k/(cp*rho)

   nx = +document.getElementById("nx").value;
   ny = +document.getElementById("ny").value;
   eps = +document.getElementById("eps").value;
   maxIt = +document.getElementById("maxIt").value;
   BC1 = +document.getElementById("BC1").value;
   BC2 = +document.getElementById("BC2").value;
   BC3 = +document.getElementById("BC3").value;
   BC4 = +document.getElementById("BC4").value;
   BC1type = String(document.getElementById("BC1type").value);
   BC2type = String(document.getElementById("BC2type").value);
   BC3type = String(document.getElementById("BC3type").value);
   BC4type = String(document.getElementById("BC4type").value);

   dx = Lx/nx;
   dy = Ly/ny;

   let T = new Array(nx).fill(0).map(()=>Array(ny).fill(0));
   let T0 = new Array(nx).fill(0).map(()=>Array(ny).fill(0));

   T = boundaryCondition1(T, k, nx, ny, dx, dy, BC1, BC1type);
   T = boundaryCondition2(T, k, nx, ny, dx, dy, BC2, BC2type);
   T = boundaryCondition3(T, k, nx, ny, dx, dy, BC3, BC3type);
   T = boundaryCondition4(T, k, nx, ny, dx, dy, BC4, BC4type);
   T = edgesBC(T, nx, ny);
   console.log(T)

   T0 = boundaryCondition1(T0, k, nx, ny, dx, dy, BC1, BC1type);
   T0 = boundaryCondition2(T0, k, nx, ny, dx, dy, BC2, BC2type);
   T0 = boundaryCondition3(T0, k, nx, ny, dx, dy, BC3, BC3type);
   T0 = boundaryCondition4(T0, k, nx, ny, dx, dy, BC4, BC4type);
   T0 = edgesBC(T0, nx, ny);

   T = solve(T, T0, nx, ny, dx, dy, BC1, BC2, BC3, BC4, BC1type, BC2type, BC3type, BC4type, eps, maxIt);
   plot(T);
}

document.getElementById("solve").addEventListener("click", compute);
  