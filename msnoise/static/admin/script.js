(function() {
	var g = {}; // graphs
	var u = {}; // utility functions
	u.padt = function(s) {
		return ('00' + s).slice(-2);
	};
	u.seconds = function(n) {
		var d = Math.floor(n / 86400);
		var h = Math.floor(n / 3600) % 24;
		var m = Math.floor(n / 60) % 60;
		var s = Math.floor(n) % 60;
		return [d, u.padt(h), u.padt(m), u.padt(s)].join(':');
	};
	u.bytes = function(n, decimal) {
		var base = decimal ? 10 : 2;
		var exp = decimal ? 3 : 10;
		var units = decimal ? ['B', 'KB', 'MB', 'GB', 'TB', 'PB'] :
		                      ['B', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB'];
		if (n < 0) {
			n = -n;
			s = '-';
		} else {
			s = '';
		}
		for (i = 5; i >= 0; i--)
			if (n >= Math.pow(base, i * exp) - 1)
				return s + (n / Math.pow(base, i * exp)).
					toFixed(2) + ' ' + units[i];
	};
	u.percent = function(n) {
		return n.toFixed(1) + '%';
	};
	var h = {}; // data handlers
	h.fqdn = function(s) {
		$('#fqdn').text(s);
	};
	h.uptime = function(n) {
		$('#uptime').text(u.seconds(n));
	};
	h.cpuusage = function(n) {
        for (var i = 0; i < n.length; i++) {
            $('#cpuusage'+i).text(u.percent(n[i]));
            g['cpuusage'.concat(i)].t.append(+new Date, n[i]);
        }
	};
	h.ramusage = function(n) {
		$('#ramusage').text(u.bytes(n[0] - n[1]));
		g.ramusage.t.append(+new Date, n[2]);
	};
	h.diskio = function(n) {
		$('#diskr').text(u.bytes(n[2], 1));
		$('#diskw').text(u.bytes(n[3], 1));
		if (h.diskio.lastr != undefined) {
			var rs = n[2] - (h.diskio.lastr || 0);
			var ws = n[3] - (h.diskio.lastw || 0);
			$('#diskrs').text(u.bytes(rs, 1) + '/s');
			$('#diskws').text(u.bytes(ws, 1) + '/s');
			g.diskrs.t.append(+new Date, rs / 1048576);
			g.diskws.t.append(+new Date, ws / 1048576);
		}
		h.diskio.lastr = n[2];
		h.diskio.lastw = n[3];
	};
	h.diskusage = function(n) {
		$('#disku').text(u.bytes(n[0], 1));
		$('#diskt').text(u.bytes(n[1], 1));
	};
	h.netio = function(n) {
		$('#nett').text(u.bytes(n[0], 1));
		$('#netr').text(u.bytes(n[1], 1));
		if (h.netio.lastt != undefined) {
			var ts = n[0] - (h.netio.lastt || 0);
			var rs = n[1] - (h.netio.lastr || 0);
			$('#netts').text(u.bytes(ts, 1) + '/s');
			$('#netrs').text(u.bytes(rs, 1) + '/s');
			g.netts.t.append(+new Date, ts / 1048576);
			g.netrs.t.append(+new Date, rs / 1048576);
		}
		h.netio.lastt = n[0];
		h.netio.lastr = n[1];
	};
	h.swapusage = function(n) {
		$('#swapusage').text(u.bytes(n[1]));
		g.swapusage.t.append(+new Date, n[3]);
	};
	var count = 0, errors = 0;
	var latency = 0;
	var wait = 1000;
	var margin = 250;
	function ping() {
		var time = +new Date;
		heartbeaton();
		$.get('/admin/rawps', function(data) {
			document.title = data.fqdn;
			for (var i in data)
				h[i] && h[i](data[i]);
			// compensate for request time while allowing time for
			// the heartbeat transitions to complete
			var t = latency = new Date - time;
			++count;
			update();
			if (t <= margin) {
				setTimeout(heartbeatoff, margin - t);
				setTimeout(ping, wait - t);
			} else if (t <= wait - margin) {
				heartbeatoff();
				setTimeout(ping, wait - t);
			} else {
				heartbeatoff();
				setTimeout(ping, margin);
			}
		});
	}
	function update() {
		$('#latency').text(latency + ' ms');
		g.latency.t.append(+new Date, latency);
		$('#requests').text(count + '/' + errors);
	}
	function heartbeaton() {
		$('#heartbeat').addClass('on');
	}
	function heartbeatoff() {
		$('#heartbeat').removeClass('on');
	}
	function error() {
		++errors;
		$('#uptime').text('OFFLINE');
		update();
		setTimeout(ping, wait);
	}
	function graph(name, percentage) {
		var options = percentage ? {
			millisPerPixel: 100,
			grid: {
				fillStyle: 'rgba(32, 32, 32, 1)',
				strokeStyle: 'rgba(0, 0, 0, 0)'
			},
			minValue: 0,
			maxValue: 100,
			labels: { fillStyle: 'rgba(0, 0, 0, 0)' }
		} : {
			millisPerPixel: 100,
			grid: {
				fillStyle: 'rgba(32, 32, 32, 1)',
				strokeStyle: 'rgba(0, 0, 0, 0)'
			}
		};
		var ts_options = {
			strokeStyle: 'rgba(255, 64, 64, 1)'
		};
		g[name] = {
			c: new SmoothieChart(options)
		};
		g[name].c.streamTo($('#c_' + name)[0], 1000);
		g[name].t = new TimeSeries();
		g[name].c.addTimeSeries(g[name].t, ts_options);
	}
	$.ajaxSetup({
		timeout: 5000,
		error: error
	});
	var graphlist = {
		latency: 0,
		ramusage: 1,
		diskrs: 0,
		diskws: 0,
		netts: 0,
		netrs: 0,
		swapusage: 1
	};
	for (var i in graphlist)
		graph(i, graphlist[i]);

    for (i=0;  i < cpucount; i++)
        graph("cpuusage".concat(i),1);

	ping();
})();
