PNGNQ - PNG NEUQUANT
===============================================

# install

    pip install pypngnq

#use

    from pypngnq import PngNQ
    n = PngNQ('/tmp/test.png')
    n.use_floyd = False
    n.exclusion_threshold = 0.3
    n.save('/tmp/test_pypngnq.png')
    n.load('/tmp/test_pypngnq.png')
    n.save('/tmp/test_pypngnq2.png')
    n.close()
