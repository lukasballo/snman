import osm2gmns as og

net = og.getNetFromFile(
    'C:\DATA\CLOUD STORAGE\polybox\Research\SNMan\SNMan Shared\inputs\osm_test_data\hdb.osm',
    POI=True,
    network_types=('auto','walk','bike')
)

og.combineShortLinks(net)

og.consolidateComplexIntersections(
    net,
    auto_identify=True,
    int_buffer=200
)

og.generateMovements(net)

og.outputNetToCSV(
    net,
    output_folder='C:\DATA\CLOUD STORAGE\polybox\Research\SNMan\SNMan Shared\outputs\osm2gmns'
)

print("Done")