
import sys

from lxml import etree

xml = etree.parse(sys.argv[1])

print """<phillip>
<configure>
<version>DUMMY</version>
<components lhs="DUM" ilp="DUM" sol="DUM"></components>
<knowledge_base path="DUM" size="-1" max_distance="-1"></knowledge_base>
<params timeout_lhs="-1" timeout_ilp="-1" timeout_sol="-1" timeout_all="-1" verbose="-1" path_out="DUM" max_depth="-1"></params>
</configure>                                                                                                                                               
"""

if len(sys.argv) == 4:
	print etree.tostring(
		xml.xpath(
			"/phillip-learn/round[@iteration='%s']/%s/proofgraph" % (
				sys.argv[2],
				sys.argv[3]))[0])

else:
	try:
		print etree.tostring(
			xml.xpath(
				"/phillip-learn/round[@iteration='%s' and @observation='%s']/%s/proofgraph" % (
					sys.argv[2],
					sys.argv[3],
					sys.argv[4]))[0])

	except IndexError:
		print ""
		print >>sys.stderr, "Oh..."
	
print "</phillip>"
