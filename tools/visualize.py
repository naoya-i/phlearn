from lxml import etree

import StringIO
import sys
import re
import collections

# USAGE: <filename> <observaiton> <itr> <fvdict>

def _createTree(xml):
	ret = []
	nodes = collections.defaultdict(list)
	n2repr   = {}
	n2prop   = {}

	def _getNodeIDs(_x):
		return _x.split(":")[-1][1:-1].split(",")

	def _deco(_x, mode):
		return "<code>%s</code>" % _x if mode else "%s" % _x

	for expl in xml.xpath("./literals/literal"):
		nodes[expl.attrib["id"]] = []
		n2repr[expl.attrib["id"]] = expl.text
		n2prop[expl.attrib["id"]] = expl

	for expl in xml.xpath("./explanations/explanation"):
		for tailNode in _getNodeIDs(expl.attrib["tail"]):
			for headNode in _getNodeIDs(expl.attrib["head"]):
				nodes[headNode] += ["e%s" % expl.attrib["id"]]

			nodes["e%s" % expl.attrib["id"]] += [tailNode]

	for node, nodeParents in nodes.iteritems():
		if 0 == len(nodeParents):
			ret += ["{id: 'n%s', text: '%s', parent: '#'}" % (
				node,
				_deco(n2repr.get(node, node), n2prop.has_key(node) and "yes" == n2prop[node].attrib["active"]),
				)]

		else:
			for np in nodeParents:
				ret += ["{id: 'n%s', text: '%s', state: {%s}, parent: 'n%s'}" % (
					node,
					_deco(n2repr.get(node, node), n2prop.has_key(node) and "yes" == n2prop[node].attrib["active"]),
					"selected: true" if n2prop.has_key(node) and "yes" == n2prop[node].attrib["active"] else "",
					np
				)]

	return ret

def html(fn, obs_name, itr, fvdict):
	m = re.search("<round iteration=\"%s\" observation=\"%s\">(.*?)</round>" % (itr, obs_name), open(fn).read(), re.DOTALL)
	xml = etree.parse(StringIO.StringIO("<phillip-learn>" + m.group(0) + "</phillip-learn>"))
	if "-" != fvdict:
		fvdict = dict([tuple(x.split("\t")) for x in open(fvdict)])
	else:
		fvdict = {}

	# Counting the list of observations and iteration.
	list_itr = xml.xpath("/phillip-learn/round[@observation='%s']/@iteration" % obs_name)
	obs      = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/current-prediction/proofgraph/literals/literal[@type='observable']/text()" % (obs_name, itr))
	labels   = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/label/text()" % (obs_name, itr))
	label_status   = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/label-status/text()" % (obs_name, itr))
	lf_cp    = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/current-prediction/logical-form/text()" % (obs_name, itr))
	lf_lvc   = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/latent-variable-completion/logical-form/text()" % (obs_name, itr))
	vec_cp    = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/current-prediction/vector/text()" % (obs_name, itr))
	vec_lvc   = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/latent-variable-completion/vector/text()" % (obs_name, itr))
	sc_cp     = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/current-prediction/proofgraph" % (obs_name, itr))[0]
	sc_lvc    = xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/latent-variable-completion/proofgraph" % (obs_name, itr))[0]

	if 0 == len(vec_cp): vec_cp = ["0:0"]
	if 0 == len(vec_lvc): vec_lvc = ["0:0"]

	sc_cp     = "%s (%s)" % (sc_cp.attrib["objective"], sc_cp.attrib["state"])
	sc_lvc    = "%s (%s)" % (sc_lvc.attrib["objective"], sc_lvc.attrib["state"])

	tree_cp = _createTree(xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/current-prediction/proofgraph" % (obs_name, itr))[0])
	tree_lvc = _createTree(xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/latent-variable-completion/proofgraph" % (obs_name, itr))[0])

	#
	lms = []

	def _node2repr(_x):
		return _x.attrib["target"] + ("" if "0" == _x.attrib["is_transitive_eq"] else " (TRANSEQ)")

	for lblno, find_label in enumerate(xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/latent-variable-completion/find-label" % (obs_name, itr))):
		lms += ["<pre>%s. %s</pre>" % (1+lblno, find_label.attrib["label"])]
		lms += ["<table class=\"table table-bordered\">"]

		lms += ["<tr>"]
		lms += ["<th>#</th>"]
		lms += ["<th>Matched nodes</th>"]
		lms += ["</tr>"]

		for patno, pat in enumerate(find_label.xpath("./pattern")):
			lms += ["<tr>"]
			lms += ["<th>%s</th>" % (1+patno)]
			lms += ["<td><pre>%s</pre></td>" % ", ".join([_node2repr(x) for x in pat.xpath("./proofgraph-checking/node")])]
			lms += ["</tr>"]

		lms += ["</table>"]

	# Weight updates.
	wupdates = {}

	for upd in xml.xpath("/phillip-learn/round[@observation='%s' and @iteration=%s]/weight-update/update" % (obs_name, itr)):
		wupdates[upd.attrib["id"]] = (upd.attrib["before_m"], upd.attrib["after_m"])

	lf_cp   = [x for x in lf_cp[0].strip().split(",") if not x.startswith("(!=")]
	lf_lvc  = [x for x in lf_lvc[0].strip().split(",") if not x.startswith("(!=")]
	vec_cp  = sorted(vec_cp[0].strip().split(" "), key=lambda x: int(x.split(":")[0]))
	vec_lvc = sorted(vec_lvc[0].strip().split(" "), key=lambda x: int(x.split(":")[0]))

	vekeys_cp = [x.split(":")[0] for x in vec_cp]
	vekeys_lvc = [x.split(":")[0] for x in vec_lvc]
	vedict_cp = dict([x.split(":") for x in vec_cp])
	vedict_lvc = dict([x.split(":") for x in vec_lvc])
	fks       = set(vekeys_cp)|set(vekeys_lvc)

	# Highlight gold features.
	def _labeling(isCorrect):
		return "<span class=\"label label-danger\">O</span>" if isCorrect else \
			"<span class=\"label label-info\">X</span>"

	def _update(x):
		if "" == x: return ""

		return "%s => %s" % x + (
			"<span class=\"label label-danger\">UP!</span>" if float(x[0]) < float(x[1]) else \
			"<span class=\"label label-info\">DOWN!</span>")


	wupdates = [
		"<tr><th>%s</th><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
			fvdict.get(fk, fk), _labeling(vedict_cp.has_key(fk) and vedict_lvc.has_key(fk)), vedict_cp.get(fk, ""), vedict_lvc.get(fk, ""), _update(wupdates.get(fk, "")),
		) for fk in fks
		]

	dct = {
		"lf_cp":        ", ".join(lf_cp),
		"lf_lvc":       ", ".join(lf_lvc),
		"vec_cp":        " ".join(vec_cp),
		"vec_lvc":       " ".join(vec_lvc),
		"obs_name":     obs_name,
		"itr":          itr,
		"obs":      ", ".join(["%s" % x for x in obs]),
		"list_itr": "".join(["<option>%s</option>" % x for x in list_itr]),
		"labels":   "".join(["%d. %s<br />" % (1+i, x) for i, x in enumerate(labels)]),
		"label_status": label_status[0],
		"tree_cp":      ",".join(tree_cp),
		"tree_lvc":      ",".join(tree_lvc),
		"weightupdate":  "".join(wupdates),
		"sc_cp":         sc_cp,
		"sc_lvc":         sc_lvc,
		"lms":           "".join(lms),
		}

	return """
<html>
	<head>
	  <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
  	<link href="./lib/bootstrap-3.0.3/dist/css/bootstrap.min.css" rel="stylesheet" />
	  <link href="./lib/vakata-jstree-a0767ce/dist/themes/default/style.min.css" rel="stylesheet" />
	  <script>
	    var g_treeCp = { 'core' : { 'data' : [%(tree_cp)s] } };
	    var g_treeLvc = { 'core' : { 'data' : [%(tree_lvc)s] } };
	  </script>
	</head>
	<body>
	<br />
	<br />
	<br />
	<div class="container-fluid" style="margin: 10px">

	<h2><span class="label label-default">Observation</span></h2>
	<pre>%(obs)s</pre>

	<h2><span class="label label-default">Labels (status: %(label_status)s)</span></h2>
	<pre>%(labels)s</pre>

	<h2><span class="label label-default">Inference result</span></h2>
	<table class="table table-bordered">
	<thead>
	<tr>
	  <th></th>
	  <th><span class="label label-success">Current Prediction</span></th>
	  <th><span class="label label-danger">Latent Variable Comp.</span></th></tr>
	</thead>
	<tbody>
	<tr><th scope="row">Logical form</th><td><pre>%(lf_cp)s</pre></td><td><pre>%(lf_lvc)s</pre></td></tr>
	<tr><th scope="row">Score</th><td>%(sc_cp)s</td><td>%(sc_lvc)s</td></tr>
	<tr><th scope="row">Feature vector</th><td>%(vec_cp)s</td><td>%(vec_lvc)s</td></tr>
	<tr><th scope="row">Proof graph</th><td><div id="jstree_cp"></div></td><td><div id="jstree_lvc"></div></td></tr>
	</tbody>
	</table>

	<h2><span class="label label-default">Weight update</span></h2>
	<table class="table table-bordered">
	<thead>
	<tr>
	  <th>#</th>
	  <th>X/O</th>
	  <th><span class="label label-success">Current Prediction</span></th>
	  <th><span class="label label-danger">Latent Variable Comp.</span></th>
	  <th>Update</th>
	</tr>
	</thead>
	<tbody>
	%(weightupdate)s
	</tbody>
	</table>


	<h2><span class="label label-default">Label matching log</span></h2>
	%(lms)s

	</div>

	<nav class="navbar navbar-inverse navbar-fixed-top">
    <div class="container">
	    <div class="navbar-header">
	      <a class="navbar-brand" href="#">ilp_loglinear.cpp::Log Tracer - %(obs_name)s/iteration %(itr)s</a>
    	</div>
	  </div>
	</nav>

	<script src="./lib/vakata-jstree-a0767ce/dist/libs/jquery.js"></script>
	<script src="./lib/vakata-jstree-a0767ce/dist/jstree.min.js"></script>
	<script src="./lib/bootstrap-3.0.3/dist/js/bootstrap.min.js"></script>

	<script>
	$(function () { $('#jstree_cp').jstree(g_treeCp); });
	$(function () { $('#jstree_lvc').jstree(g_treeLvc); });
	</script>

	</body>
</html>""" % dct

if "__main__" == __name__:
	print html(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
