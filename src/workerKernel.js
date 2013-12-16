/**
 * @author Petar Petrov
 * Communication via channels: http://threadcomm.appspot.com/
 
 */
 
 /*
var channelPort;
onmessage = function(e) {
	if (e.data.code == "start") {
		channelPort = e.ports[0];
		channelPort.postMessage(e.data.msg+" >> worker1 channel post");
		channelPort.onmessage = getChannelMessage;
	} else if(e.data.code=="msgw")
	{
		postMessage(e.data.msg+" >> worker1 got msg");
	}else if(e.data.code=="msgch")
	{
		channelPort.postMessage(e.data.msg+" >> worker1 got msg");
	}
}
function getChannelMessage(e){
	postMessage(e.data+" >> channel recieved msg in worker1 ");
}
*/
self.onmessage = function(e) 
{
	log("Web Worker :: START");
	
	if (e.data.url) {
		var url = e.data.url;
		var urlindex = url.indexOf('index.html');
		if (urlindex != -1) {
			url = url.substring(0, urlindex);
		}
		importScripts(url+"libs/numeric/numeric-1.2.6.js");
	}
	
	if (e.data.code == "sim") 
	{
		var params = e.data.params;
		var numPointsDomain = e.data.points.length;
		var numPointsCurv = e.data.curve.length;
		var L = e.data.curve;
		var results = new Array(numPointsDomain);
		
		//for each point in the input calculate induction
		for(i = 0; i < numPointsDomain; i++)
		{
			var point = e.data.points[i];
			results[i] = [0,0,0];
			
			for(c = 0; c < numPointsCurv-1; c++)
			{
				//Length of the curve element
				var diflen = numeric.sub(L[c],L[c+1]);
				var len = numeric.norm2(diflen);

				//Number of points for the curve-element discretization
				var len_ds = numeric.div(len,params.ds);
				var Npi = Math.ceil(len/params.ds);
				if(Npi < 3){
					log("ERROR Integration step is too big!!");
				}

				//% Curve-element discretization
				var Lx = numeric.linspace(L[c][0], L[c+1][0], Npi);
				var Ly = numeric.linspace(L[c][1], L[c+1][1], Npi);
				var Lz = numeric.linspace(L[c][2], L[c+1][2], Npi);

				var Ldiscrete = new Array(Npi);
				for(ci = 0;ci< Npi; ci++)
				{
					Ldiscrete[ci] = [Lx[ci],Ly[ci],Lz[ci]];
				}

				//% Integration
				for(s = 0;s<Npi-1;s++)
				{
					//Vector connecting the infinitesimal curve-element			
					var Rxyz = numeric.sub(Ldiscrete[s] , point);

					//Infinitesimal curve-element components
					var dLxyz = numeric.sub(Ldiscrete[s+1] , Ldiscrete[s]);

					//Modules
					var dLn = numeric.norm2(dLxyz);
					var Rn = numeric.norm2(Rxyz);

					//Biot-Savart
					var dB = numeric.mul( crossprod(Rxyz,dLxyz), 
						params.gamma/4/params.pi/Rn/Rn/Rn);

					//Debug
					//log("BS dB: "+dB);
						
					//Add increment to the main field
					results[i] = numeric.add(results[i], numeric.mul(dB, params.loops) );
				}
			}
			
			//Debug
			//log("result [i]:value ::"+i+":"+results[i]);
			
			self.postMessage({code:'progress'});
		}
		
		self.postMessage({ID: e.data.ID, code:'sim', output:results});
		self.postMessage({ID: e.data.ID, code:'finished'});
	}
	
	log("Web Worker :: FINISH");
}

//should be in numeric.js already, but NOT :(
function crossprod(u,v)
{
	return [u[1]*v[2]-u[2]*v[1],
			u[2]*v[0]-u[0]*v[2],
			u[0]*v[1]-u[1]*v[0]];
}

function log(msg)
{
	self.postMessage({code:'debug',message:msg});
}