//Author: Petar Petrov
//Credits: Alessandro Masullo
//Matlab Source: http://www.mathworks.nl/matlabcentral/fileexchange/42237-biot-savart-direct-integration-on-a-generic-curve/content/biot_savart.m
//

var ProtoMag = namespace("DeNaCoR.demo1");

ProtoMag.Sim = function(args)
{
	//draw grid
	//ProtoMag.drawGrid();

	//var progress = 1;			
	ProtoMag.progress(1);

	//define domain
	var ND = args.dim;
	//var Dom = 	[[-1.1, 1],
	//			[-1.1, 1],
	//			[0.1, 10]];
	
	var Dom = 	[[-args.stepx/2 * args.dim, args.stepx/2 * args.dim],
				[-args.stepy/2 * args.dim, args.stepy/2 * args.dim],
				[-args.stepz/2 * args.dim, args.stepz/2 * args.dim]];
	
	// Induction constant
	var gamma = 1;

	//% Integration step size
	var ds = 0.1;

	//%% Induction curve
	var coil = ProtoMag.circleGen(0,0,0,0.5, args.segs);
	var L;
	var Nl;
	var menu1selection = 4;
	switch(menu1selection) //menu('Choose a test case:', 'Straight wire', 'Bent wire', 'Solenoid');
	{
	case 1:
		//% Test case 1: Straight wire
		L = [[0, 0, 0],
			 [0, 0, 5]];
		Nl = 2;//sizof(L)/3
		break;
	case 2:
		//% Test case 2: Bent wire
		L = [[-.5, 0, 0],
			 [-.5, 0, 4],
			 [1.5, 0, 4]];
		Nl = 3;//sizeof(L)/3
		break;
	case 3:
		//% Test case 3: Solenoid
		var theta = numeric.linspace(0, 15*pi, 70);
		//L = [cos(theta') sin(theta') theta'/10];
		break;
	case 4:
		L = coil;
		Nl = args.segs;
		break;
	}

	//%% Declaration of variables
	//% Induction vector components B = (U, V, W);
	var UVW = ProtoMag.createArray(ND, ND, ND, [0,0,0]);

	//% Volume Mesh
	var mesh = ProtoMag.meshgrid(
						numeric.linspace(Dom[0][0], Dom[0][1], ND),
						numeric.linspace(Dom[1][0], Dom[1][1], ND),
						numeric.linspace(Dom[2][0], Dom[2][1], ND));
	var step = 1;
	var final_step = (ND*ND*ND);
	
	
	for(i = 0; i<ND; i++){
		for(j = 0;j<ND;j++){
			for(k = 0;k<ND;k++){

			//% Ptest is the point of the field where we calculate induction
			var pTest = mesh[i][j][k];//mesh[X(i,j,k),Y(i,j,k),Z(i,j,k)];

			//% The curve is discretized in Nl points, we iterate on the Nl-1
			//% segments. Each segment is discretized with a "ds" length step
			//% to evaluate a "dB" increment of the induction "B".
			for(pCurv = 0;pCurv<Nl-1;pCurv++)
			{
				//% Length of the curve element
				var diflen = numeric.sub(L[pCurv],L[pCurv+1]);
				var len = numeric.norm2(diflen);

				//% Number of points for the curve-element discretization
				var len_ds = numeric.div(len,ds);
				var Npi = Math.ceil(len/ds);
				if(Npi < 3){
					//ProtoMag.log("ERROR Integration step is too big!!");
				}

				//% Curve-element discretization
				var Lx = numeric.linspace(L[pCurv][0], L[pCurv+1][0], Npi);
				var Ly = numeric.linspace(L[pCurv][1], L[pCurv+1][1], Npi);
				var Lz = numeric.linspace(L[pCurv][2], L[pCurv+1][2], Npi);

				var Ldiscrete = new Array(Npi);
				for(ii = 0;ii< Npi; ii++)
				{
					Ldiscrete[ii] = [Lx[ii],Ly[ii],Lz[ii]];
				}

				//% Integration
				for(s = 0;s<Npi-1;s++)
				{
					//% Vector connecting the infinitesimal curve-element
					//% point and field point "pTest"
					var Rxyz = numeric.sub(Ldiscrete[s] , pTest);
					//var Rx = numeric.sub(Lx[s] , pTest[0]);
					//var Ry = numeric.sub(Ly[s] , pTest[1]);
					//var Rz = numeric.sub(Lz[s] , pTest[2]);

					//% Infinitesimal curve-element components
					var dLxyz = numeric.sub(Ldiscrete[s+1] , Ldiscrete[s]);
					//var dLx = numeric.sub(Lx[s+1] , Lx[s]);
					//var dLy = numeric.sub(Ly[s+1] , Ly[s]);
					//var dLz = numeric.sub(Lz[s+1] , Lz[s]);

					//% Modules
					//var dL = Math.sqrt(Math.pow(dLx,2) + Math.pow(dLy,2) + Math.pow(dLz,2));
					var dLn = numeric.norm2(dLxyz);
					//var R = Math.sqrt(Math.pow(Rx,2) + Math.pow(Ry,2) + Math.pow(Rz,2));
					var Rn = numeric.norm2(Rxyz);

					//% Biot-Savart
					var dUVW = numeric.mul( ProtoMag.crossprod(Rxyz,dLxyz), gamma/4/pi/Rn/Rn/Rn);
					//var dU = gamma/4/pi*(dLy*Rz - dLz*Ry)/R/R/R;
					//var dV = gamma/4/pi*(dLz*Rx - dLx*Rz)/R/R/R;
					//var dW = gamma/4/pi*(dLx*Ry - dLy*Rx)/R/R/R;

					//% Add increment to the main field
					UVW[i][j][k] = numeric.add(UVW[i][j][k], numeric.mul(dUVW, args.loops) );
					//U[i][j][k] += dU;
					//V[i][j][k] += dV;
					//W[i][j][k] += dW;
				}
			}

			ProtoMag.progress( Math.floor(  ( step / final_step ) * 100 ));
			step++;
		}}}

		ProtoMag.progress(100);


		ProtoMag.CoildMesh = coil;
		ProtoMag.Mesh = mesh;
		ProtoMag.VectorField = UVW;
}

ProtoMag.showVectorField = function(node,fieldType,normalized,vecPos,vecDir,cutZ)
{
	var xgvL = vecPos.length;
	var ygvL = vecPos[0].length;
	var zgvL = vecPos[0][0].length;
	
	var vecNorm = ProtoMag.normalizeArray3(vecDir);
	
	var MinMax = ProtoMag.getMinMaxField(vecDir);			
	

	
	for(i = 0; i < xgvL; i++)
	{
		for(j = 0; j < ygvL; j++)
		{
			for(k = 0; k < zgvL; k++)
			{
				if(cutZ < 1 || k == cutZ)
				{
				
				
					var mesh = null;
													
					var mag = numeric.norm2(vecDir[i][j][k]);
					
					var color = ProtoMag.getColorMapValue(MinMax,mag);
					
					if(normalized)
						mag = 0.25;
					else
						mag = numeric.norm2(vecDir[i][j][k]);
						
					
					
					
					if(fieldType == 0)
					{
					   //( radiusTop, radiusBottom, height, radialSegments, heightSegments, openEnded )
						var geometry = new THREE.CylinderGeometry(0,.05,0.4,4,2,false);
						var material = new THREE.MeshBasicMaterial( { color: color } );
						
						mesh = new THREE.Mesh( geometry, material );

						mesh.position = new THREE.Vector3(vecPos[i][j][k][0],
														vecPos[i][j][k][1],
														vecPos[i][j][k][2]);
														
						ProtoMag.setDir(mesh, 
							new THREE.Vector3(vecNorm[i][j][k][0],
											vecNorm[i][j][k][1],
											vecNorm[i][j][k][2]));
						
					}			
					else if(fieldType == 1)
					{
						mesh = new THREE.ArrowHelper(
							new THREE.Vector3(vecNorm[i][j][k][0],vecNorm[i][j][k][1],vecNorm[i][j][k][2]), 
							new THREE.Vector3(vecPos[i][j][k][0],vecPos[i][j][k][1],vecPos[i][j][k][2]),
							mag + 0.01,
							color);
					}				
					else if(fieldType == 2)
					{	

						var rad = ( mag / MinMax.max ) * 0.1 + 0.01 ;// Math.min(mag * 10.5, 1.0);
						
						//radius segments ring
						var geometry = new THREE.SphereGeometry( rad , 8, 8);
						var material = new THREE.MeshBasicMaterial( { color: color } );
						
						//centre, radius
						var mesh = new THREE.Mesh( geometry, material);
						mesh.position = new THREE.Vector3(vecPos[i][j][k][0],vecPos[i][j][k][1],vecPos[i][j][k][2]);
					}
					
					node.add( mesh );
				}
			}
		}
	}
	
	
}


ProtoMag.getColorMapValue = function (MinMax, val)
{
			
	if(ProtoMag.ColorMap === undefined)
	{
		ProtoMag.ColorMap = [
		//Blues
		0x0000FF,
		0x005CFF,
		0x00B9FF,
		0x00FFE7,
		//Greens
		0x00FF8B,
		0x00FF2E,
		0x2EFF00,
		0x8BFF00,
		//Yellow-Orange-Red
		0xE7FF00,
		0xFFB900,
		0xFF5C00,
		0xFF0000
		];
	}

	var range = MinMax.max - MinMax.min;
	
	if( range <= 0)
		ProtoMag.log("ERROR minmax mismatch");
	
	var step = range / ProtoMag.ColorMap.length;
	
	var index =  val / step;
	
	index = Math.min( Math.floor(index), ProtoMag.ColorMap.length - 1 );
	
	//DEBUG
	if(index < 0 || index > 11 )
		ProtoMag.log("WRONG ColorMap index: "+index);
	
	return ProtoMag.ColorMap[index];			
}
   
ProtoMag.setDir = function (mesh,dir)
{
var axis = new THREE.Vector3();
var radians;
// dir is assumed to be normalized

	if ( dir.y > 0.99999 ) {

			mesh.quaternion.set( 0, 0, 0, 1 );

	} else if ( dir.y < - 0.99999 ) {

			mesh.quaternion.set( 1, 0, 0, 0 );

	} else {

			axis.set( dir.z, 0, - dir.x ).normalize();

			radians = Math.acos( dir.y );

			mesh.quaternion.setFromAxisAngle( axis, radians );
	}
}

ProtoMag.showLineSegments = function(node,pointwise)
{
	var geometry = new THREE.Geometry();
	var material = new THREE.LineBasicMaterial( { color: 0xEEEEEE, linewidth : 3 } );

	for(i = 0; i < pointwise.length; i++)
	{
		geometry.vertices.push(new THREE.Vector3(pointwise[i][0],pointwise[i][1],pointwise[i][2]));
	}

	var lines = new THREE.Line( geometry, material);
	lines.position.set(0,0,0);

	node.add(lines);
	
	//LineBasicMaterial( parameters )
	//parameters is an object with one or more properties defining the material's appearance.
	//color — Line color in hexadecimal. Default is 0xffffff.
	//linewidth — Line thickness. Default is 1.
	//linecap — Define appearance of line ends. Default is 'round'.
	//linejoin — Define appearance of line joints. Default is 'round'.
	//vertexColors — Define whether the material uses vertex colors, or not. Default is false.
	//fog — Define whether the material color is affected by global fog settings. Default is false.
	//%% Graphic
	//figure(1)
}

ProtoMag.drawGrid = function()
{
	var canvas = document.getElementById("canvas");
	var context = canvas.getContext("2d");
	//b_context.fillRect(50, 25, 150, 100);

	console.log('canvas.height:'+canvas.height);

	for (var x = 0.5; x < canvas.width; x += 10) {
	  context.moveTo(x, 0);
	  context.lineTo(x, canvas.height);
	}

	for (var y = 0.5; y < canvas.height; y += 10) {
	  context.moveTo(0, y);
	  context.lineTo(canvas.width, y);
	}
	context.strokeStyle = "#eee";
	context.stroke();

}


ProtoMag.meshgrid = function(xgv,ygv,zgv)
{
	//[X,Y,Z] = meshgrid(xgv,ygv,zgv) produces three-dimensional coordinate arrays.
	//The output coordinate arrays X, Y, and Z contain copies of the grid vectors xgv, ygv, and zgv respectively.
	//The sizes of the output arrays are determined by the length of the grid vectors.
	//For grid vectors xgv, ygv, and zgv of length M, N, and P respectively, X, Y, and Z will have N rows, M columns, and P pages.

   var xgvL = xgv.length;
   var ygvL = ygv.length;
   var zgvL = zgv.length;

   var M = ProtoMag.createArray(xgvL, ygvL, zgvL, [0,0,0]);

   for(i = 0; i < xgvL; i++)
   {
	   for(j = 0; j < ygvL; j++)
	   {
			for(k = 0; k < zgvL; k++)
			{
				M[i][j][k] = [xgv[i],ygv[j],zgv[k]];
			}
	   }
   }

	return M;
}

ProtoMag.createArray = function(a,b,c,val)
{
   var A = new Array(a);

   for(i = 0; i<a; i++)
   {
	   A[i] = new Array(b);

	   for(j = 0; j < b; j++)
	   {
		   A[i][j] = new Array(c);

		   for(k = 0; k < c; k++)
		   {
			   A[i][j][k] = val;
		   }
	   }
   }
   return A;
}

ProtoMag.normalizeArray3 = function(array)
{
	var dim1 = array.length;
	var dim2 = array[0].length;
	var dim3 = array[0][0].length;

	var narray = ProtoMag.createArray(dim1,dim2,dim3,[0,0,0]);

   for(i = 0; i < dim1; i++)
   {
	   for(j = 0; j < dim2; j++)
	   {
			for(k = 0; k < dim3; k++)
			{
			 var n = numeric.norm2(array[i][j][k]);
			 narray[i][j][k] = numeric.div(array[i][j][k], n);
			}
		}
	}

	return narray;
}
   
ProtoMag.getMinMaxField = function(array)
{
	var dim1 = array.length;
	var dim2 = array[0].length;
	var dim3 = array[0][0].length;

	var max = 0;
	var min = 9999999;

   for(i = 0; i < dim1; i++)
   {
	   for(j = 0; j < dim2; j++)
	   {
			for(k = 0; k < dim3; k++)
			{
			 var val = numeric.norm2(array[i][j][k]);
			 
			 if(val > max) max = val;
			 
			 if(val < min) min = val;
			 
			}
		}
	}

	return {max : max, min : min};
}

ProtoMag.log = function(text)
{
   console.log(text);
}

ProtoMag.progress = function(percent)
{
   //$("#progressbar").progressbar({ value: percent });
   setProgress(percent);
}

ProtoMag.crossprod = function(u,v)
{
	return [u[1]*v[2]-u[2]*v[1],
			u[2]*v[0]-u[0]*v[2],
			u[0]*v[1]-u[1]*v[0]];
}

ProtoMag.circleGen = function(x,y,z,r,nsegments)
{
	//(x + r*cos(alpha), y + r*sin(alpha)
	var pi = Math.atan(1)*4;
	var anglePerSegment = 2*pi/nsegments;
	var circlexy = new Array(nsegments+1);

	for(i = 0; i < nsegments; i++)
	{
		//circlexy[i] = new THREE.Vector3(x + r * Math.cos(anglePerSegment*i), y + r * Math.sin(anglePerSegment*i), z);
		circlexy[i] = [x + r * Math.cos(anglePerSegment*i), y + r * Math.sin(anglePerSegment*i), z];
	}

	circlexy[nsegments] = circlexy[0];

	return circlexy;
}
