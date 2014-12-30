//Author: Petar Petrov
//Credits: Alessandro Masullo MATLAB Biot-Savart Implementation
//Matlab Source: http://www.mathworks.nl/matlabcentral/fileexchange/42237-biot-savart-direct-integration-on-a-generic-curve/content/biot_savart.m
//

var ProtoMag = namespace("DeNaCoR.demo1");


ProtoMag.SimWithWorkers = function(model)
{		
	this.progress(1);

	var i,j,k;
	

	//define domain
	var ND = model.domain.dimension;
	var stepx = model.domain.stepx;
	var stepy = model.domain.stepy;
	var stepz = model.domain.stepz;
	
	//var Dom = 	[[-1.1, 1],
	//			[-1.1, 1],
	//			[0.1, 10]];
	
	var Dom = 	[[-stepx/2 * ND, stepx/2 * ND],
				[-stepy/2 * ND, stepy/2 * ND],
				[-stepz/2 * ND, stepz/2 * ND]];
	
	// Induction constant
	var gamma = 1;

	//NOW this is auto
	//% Integration step size
	var ds = 0.1;

	//%% Induction curve
	var coil = ProtoMag.circleGen(0,0,0,0.5, model.coil.segments);
	var L = coil;
	var Nl = model.coil.segments;


	//%% Declaration of variables
	// Induction vector components B = (U, V, W);
	//var UVW = ProtoMag.createArray(ND, ND, ND, [0,0,0]);
	var B = new Array(ND*ND*ND);

	// Regular Grid Mesh (Sample Point Cloud for which we derive B)
	var mesh = ProtoMag.meshgrid(
						numeric.linspace(Dom[0][0], Dom[0][1], ND),
						numeric.linspace(Dom[1][0], Dom[1][1], ND),
						numeric.linspace(Dom[2][0], Dom[2][1], ND));
	var step = 1;
	var final_step = (ND*ND*ND);
	var steps_per_worker = final_step / 4;
	var nn = 0;
	var pp = 0;
	var workerIDX = 0;
	var workerData = new Array(4);
	
	var workerResults = new Array(4);
	var workers = new Array(4);
	
	for(i = 0; i<4; i++)
	{
		workerData[i] = new Array(); 
		workers[i] = new Worker("src/workerKernel.js");
	}
	
	
	//decompose input for parallel execute
	for(i = 0; i<ND; i++){
		for(j = 0;j<ND;j++){
			for(k = 0;k<ND;k++){
				
				nn++;
				
				if(nn > steps_per_worker)
				{
					workerIDX++;
					nn = 0;
				}
				
				workerData[workerIDX].push(mesh[i][j][k]);
			}
		}
	}


	for(i = 0; i<4; i++)
	{
		workers[i].postMessage({
			ID: i,
			points: workerData[i],
			curve: coil, 
			params: {
					loops: model.coil.loops,
					gamma: gamma,
					ds : ds,
					pi : pi
					},
			code: 'sim',
			url: document.location.href
			});
		
		workers[i].onmessage = function(e)
		{
			if(e.data.code == 'sim')
			{
				workerResults[e.data.ID] = e.data.output;						
			}
			
			if ( e.data.code === 'finished' ) {
                
				this.finished = true;
				
				for(i = 0; i<4; i++)
				{
					if(!workers[i].finished)
						return;
				}
				
					//compose output from parallel execute
				var nn = 0;
				var workerIDX = 0;
				for(i = 0; i<ND; i++){
					for(j = 0;j<ND;j++){
						for(k = 0;k<ND;k++){
							
							UVW[i][j][k] = workerResults[workerIDX][nn];
							
							nn++;
						
							if(nn > steps_per_worker)
							{
								workerIDX++;
								nn = 0;
							}							
						}
					}
				}
				
				ProtoMag.CoildMesh = coil;
				ProtoMag.Mesh = mesh;
				ProtoMag.VectorField = UVW;
				
				ProtoMag.Vis(model);
			}
			
			if(e.data.code == 'debug')
				console.log(e.data.message);
			
			if(e.data.code == 'progress')
			{
				//this.progress(pp++);
				var p = Math.floor( ( step++ / final_step ) * 100 );
				
				if( p != pp)
				{
					pp = p;
					ProtoMag.progress(p);
					console.log("progress: "+p);
				}
			}
			
		};
		workers[i].onerror = function(event)
		{
			throw new Error(event.message + " (" + event.filename + ":" + event.lineno + ")");
		};
	}
}

ProtoMag.Vis = function(model)
{
	var canvasElement = document.getElementById("canvas");

	this.colormap = model.colormap;
	this.scene = new THREE.Scene();

	//var camera = null;
	if(model.camera.type == 0)
		this.camera = new THREE.OrthographicCamera( -5, 5, 5, -5, 0.1, 1000 );
	else
		this.camera = new THREE.PerspectiveCamera( 75, canvasElement.clientWidth / canvasElement.clientHeight, 0.1, 100 );

	this.camera.position.x = 0;
	this.camera.position.y = 0;
	this.camera.position.z = 15;
	
	//var renderer = null;
	if (window.WebGLRenderingContext && !model.ForceDisableWebGL)
	{
		this.renderer = new THREE.WebGLRenderer({canvas:canvasElement});
	}
	else
	{
		this.renderer = new THREE.CanvasRenderer({canvas:canvasElement});		
	}
	this.renderer.setSize( canvasElement.clientWidth, canvasElement.clientHeight );
	
	//var mesh = ProtoMag.Mesh;
	//var UVW = ProtoMag.VectorField;
	//var coil = ProtoMag.CoildMesh;

	if(this.Mesh === undefined || this.VectorField === undefined)
	{
		ProtoMag.log('cannot run visualization, run simulation first');
		return;
	}

	this.rootNode = new THREE.Object3D();
	
	//show field
	ProtoMag.showVectorField(
					this.rootNode,
					model.field.type,
					model.field.normalized,
					this.Mesh,
					this.VectorField,
					model.field.cutter);
		
	//show coil
	ProtoMag.showLineSegments(this.rootNode,this.CoildMesh);
	
	//ProtoMag.root = rootNode;	
	this.scene.add(this.rootNode);
	
	
	//http://threejs.org/examples/misc_controls_trackball.html
	var controls = new THREE.TrackballControls(this.camera,this.renderer.domElement);
	controls.rotateSpeed = 1.0;
	controls.zoomSpeed = 1.2;
	controls.panSpeed = 0.8;
	controls.noZoom = false;
	controls.noPan = false;
	controls.staticMoving = true;
	controls.dynamicDampingFactor = 0.3;
	controls.keys = [ 65, 83, 68 ];
	//controls.addEventListener( 'change', Render );
	this.controls = controls;
	
	
	//Render Scene			
	//ProtoMag.renderer = renderer;
	//ProtoMag.camera = camera;
	//ProtoMag.scene = scene;			
	this.Animate();
}
ProtoMag.ShowVectorFieldCutz = function(model)
{
	this.scene.remove(this.rootNode);
	this.rootNode = new THREE.Object3D();
	//show field
	ProtoMag.showVectorField(
					this.rootNode,
					model.field.type,
					model.field.normalized,
					this.Mesh,
					this.VectorField,
					model.field.cutter);
		
	//show coil
	ProtoMag.showLineSegments(this.rootNode,this.CoildMesh);
	this.scene.add(this.rootNode);
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

/*
					if(color === undefined)
					{
						ProtoMag.log("FATAL ERROR (possible no field or wrongn one");
						return;
					}
*/

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
							new THREE.Vector3(vecNorm[i][j][k][0],
											vecNorm[i][j][k][1],
											vecNorm[i][j][k][2]), 
							new THREE.Vector3(vecPos[i][j][k][0],
											vecPos[i][j][k][1],
											vecPos[i][j][k][2]),
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
						mesh.position = new THREE.Vector3(vecPos[i][j][k][0],
														vecPos[i][j][k][1],
														vecPos[i][j][k][2]);
					}
					
					node.add( mesh );
				}
			}
		}
	}
	
	
}


ProtoMag.getColorMapValue = function (MinMax, val)
{
	
	var range = MinMax.max - MinMax.min;
	
	if( range <= 0)
		ProtoMag.log("ERROR minmax mismatch");
	
	var step = range / this.colormap.length;
	
	var index =  val / step;
	
	index = Math.min( Math.floor(index), this.colormap.length - 1 );
	
	//DEBUG
	if(index < 0 || index > 11 )
		ProtoMag.log("WRONG ColorMap index: "+index);
	
	return this.colormap[index];
}

// dir is assumed to be normalized   
ProtoMag.setDir = function (mesh,dir)
{
	var axis = new THREE.Vector3();
	var radians;

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



ProtoMag.meshgrid = function(xgv,ygv,zgv)
{
	//[X,Y,Z] = meshgrid(xgv,ygv,zgv) produces three-dimensional coordinate arrays.
	//The output coordinate arrays X, Y, and Z contain copies of the grid vectors xgv, ygv, and zgv respectively.
	//The sizes of the output arrays are determined by the length of the grid vectors.
	//For grid vectors xgv, ygv, and zgv of length M, N, and P respectively, X, Y, and Z will have N rows, M columns, and P pages.

   var xgvL = xgv.length;
   var ygvL = ygv.length;
   var zgvL = zgv.length;
   var SIZE = xgvL*ygvL*zgvL;

   //var M = ProtoMag.createArray(xgvL, ygvL, zgvL, [0,0,0]);
   

   // for(i = 0; i < xgvL; i++)
   // {
	  //  for(j = 0; j < ygvL; j++)
	  //  {
			// for(k = 0; k < zgvL; k++)
			// {
			// 	M[i][j][k] = [xgv[i],ygv[j],zgv[k]];
			// }
	  //  }
   // }

   var M = new Array(SIZE);

   for (var i = 0; i < SIZE; i++) {
   		M[i][j][k] = [xgv[i],ygv[j],zgv[k]];
   };
   
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

 ProtoMag.Animate = function()
{
	ProtoMag.Render();
	requestAnimationFrame(ProtoMag.Animate);
}
		
ProtoMag.Render = function()
{
	ProtoMag.renderer.render(ProtoMag.scene, ProtoMag.camera);
	ProtoMag.controls.update();
}
