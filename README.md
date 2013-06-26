### What is PyPngnq?

PyPngnq is a tool for quantizing PNG images in RGBA format

### Installation

You need to install libpng first, To install PyPngnq:

Debian:

    apt-get install libpng-dev -y

FreeBSD:

    cd /usr/ports/graphics/png
    make install clean

To install PyPngnq:

    pip install PyPngnq

### How to use it?

Following is a simple example:

    # -*- coding: utf8 -*-
    from pypngnq import PngNQ
    nq = PngNQ('/tmp/test.png')
    #print nq.config
    nq.use_floyd = False
    nq.save('/tmp/test_pypngnq.png')
    nq.close()
