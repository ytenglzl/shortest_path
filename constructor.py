#!/usr/bin/env python

""" Author: Yue Teng

    It gets output of rgen as its input, and construct an unweighted,
    undirected graph. The output is of the following format:

    V 5
    E {<1,2>,<2,4>,<0,4>,<3,2>}

    which means that there are 5 vertices indexed from 0 to 4, and the
    edges are presented following E with pairs of vertices. 
"""

import sys
import re
import numpy as np
import copy


def gnrtseg(x1, y1, x2, y2, name):
    """Generate line function with
    scope in x axis or y axis.
    """

    if x1==x2:
        return ((1., 0.), x1, (x1, x1), (min(y1, y2), max(y1, y2)),
            (x1, y1), (x2, y2), name)
    elif y1==y2:
        return ((0., 1.), y1, (min(x1, x2), max(x1, x2)), (y1, y1),
            (x1, y1), (x2, y2), name)
    else:
        dx = x2-x1
        dy = y2-y1
        ct = dy*x1-dx*y1
        A = (dy, -dx)
        b = ct
        return (A, b, (min(x1, x2), max(x1, x2)),
            (min(y1, y2), max(y1, y2)),
            (x1, y1), (x2, y2), name)

def calc_itsc(sgm0, sgm1):
    """Output the intersection if the
    two segments cross each other in
    the closed scope defined.
    """

    A = np.array([sgm0[0], sgm1[0]])
    b = np.array([sgm0[1], sgm1[1]])
    if np.linalg.det(A)==0: # singular
        # A = [[1, 0], [1, 0]]
        if sgm0[0][1]==0 and sgm1[0][1]==0:
            # on the same line
            if sgm0[1]==sgm1[1]:
                # no overlap
                if sgm0[3][1]<sgm1[3][0]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                elif sgm0[3][0]>sgm1[3][1]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                # overlap
                else:
                    ylims = [sgm0[3][0], sgm0[3][1],
                             sgm1[3][0], sgm1[3][1]]
                    ylims.sort()
                    return (np.array([sgm0[2][0], ylims[1]]),
                            np.array([sgm0[2][0], ylims[2]]))
            else:
                return (np.array([None, None]),
                        np.array([None, None]))
        # A = [[0, 1], [0, 1]]
        elif sgm0[0][0]==0 and sgm1[0][0]==0:
            # on the same line
            if sgm0[1]==sgm1[1]:
                # no overlap
                if sgm0[2][1]<sgm1[2][0]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                elif sgm0[2][0]>sgm1[2][1]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                # overlap
                else:
                    xlims = [sgm0[2][0], sgm0[2][1],
                             sgm1[2][0], sgm1[2][1]]
                    xlims.sort()
                    return (np.array([xlims[1], sgm0[3][0]]),
                            np.array([xlims[2], sgm0[3][0]]))
            else:
                return (np.array([None, None]),
                        np.array([None, None]))
        else:
            # on the same line
            if sgm0[1]==0 and sgm1[1]==0:
                # no overlap
                if sgm0[2][1]<sgm1[2][0]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                elif sgm0[2][0]>sgm1[2][1]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                # overlap
                else:
                    # sign of slope
                    if sgm0[0][0]/sgm0[0][1]<0:
                        xlims = [sgm0[2][0], sgm0[2][1],
                                 sgm1[2][0], sgm1[2][1]]
                        xlims.sort()
                        ylims = [sgm0[3][0], sgm0[3][1],
                                 sgm1[3][0], sgm1[3][1]]
                        ylims.sort()
                        return (np.array([xlims[1], ylims[1]]),
                                np.array([xlims[2], ylims[2]]))
                    else:
                        xlims = [sgm0[2][0], sgm0[2][1],
                                 sgm1[2][0], sgm1[2][1]]
                        xlims.sort()
                        ylims = [sgm0[3][0], sgm0[3][1],
                                 sgm1[3][0], sgm1[3][1]]
                        ylims.sort(reverse = True)
                        return (np.array([xlims[1], ylims[1]]),
                                np.array([xlims[2], ylims[2]]))
            elif (sgm0[1]!=0 and sgm1[1]!=0 and
                  (np.array(sgm0[0])/np.array(sgm1[0]))[0]==(sgm0[1]/sgm1[1])):
                # no overlap
                if sgm0[2][1]<sgm1[2][0]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                elif sgm0[2][0]>sgm1[2][1]:
                    return (np.array([None, None]),
                            np.array([None, None]))
                # overlap
                else:
                    # sign of slope
                    if sgm0[0][0]/sgm0[0][1]<0:
                        xlims = [sgm0[2][0], sgm0[2][1],
                                 sgm1[2][0], sgm1[2][1]]
                        xlims.sort()
                        ylims = [sgm0[3][0], sgm0[3][1],
                                 sgm1[3][0], sgm1[3][1]]
                        ylims.sort()
                        return (np.array([xlims[1], ylims[1]]),
                                np.array([xlims[2], ylims[2]]))
                    else:
                        xlims = [sgm0[2][0], sgm0[2][1],
                                 sgm1[2][0], sgm1[2][1]]
                        xlims.sort()
                        ylims = [sgm0[3][0], sgm0[3][1],
                                 sgm1[3][0], sgm1[3][1]]
                        ylims.sort(reverse = True)
                        return (np.array([xlims[1], ylims[1]]),
                                np.array([xlims[2], ylims[2]]))
            else:
                return (np.array([None, None]),
                        np.array([None, None]))
    else:
        x = np.linalg.solve(A, b)
        # intersection in closed scope
        error = 0.00000000000002 # round-off error
        if (x[0]>=sgm0[2][0]-error and x[0]<=sgm0[2][1]+error
            and x[1]>=sgm0[3][0]-error and x[1]<=sgm0[3][1]+error
            and x[0]>=sgm1[2][0]-error and x[0]<=sgm1[2][1]+error
            and x[1]>=sgm1[3][0]-error and x[1]<=sgm1[3][1])+error:
            return (x, np.array([None, None]))
        else:
            return (np.array([None, None]),
                    np.array([None, None]))

def gnrtvtxs(sgpr_itsc):
    """Generate all the vertexes with labels,
    with segment-pair-intersection dictionary
    as input.
    """

    itscs = [] # list storing intersections np.arrays
    ends = [] # list storing ends tuples
    for sgmpr, itsc in sgpr_itsc.items():
        # If intersection exists, vertexes exist.
        if tuple(itsc[0])!=(None, None):
            itscs.append(tuple(itsc[0]))
            if tuple(itsc[1])!=(None, None):
                itscs.append(tuple(itsc[1]))
            ends.append(sgmpr[0][4])
            ends.append(sgmpr[0][5])
            ends.append(sgmpr[1][4])
            ends.append(sgmpr[1][5])
    # discard the repeated ends and intersections
    itscs = set(itscs)
    itscs = list(itscs)
    ends = set(ends)
    ends = list(ends)

    return (itscs, ends)

def inrtitscs(sgpr_itsc, strts):
    """insert all the intersections to streets.
    """

    # a deep deplicate of strts for inserted streets
    # so that strts is not inserted
    strts_ = copy.deepcopy(strts)
    # insert
    for sgmpr, itsc in sgpr_itsc.items():
        itsc0tp = tuple(itsc[0])
        itsc0np = np.array(itsc[0])
        itsc1tp = tuple(itsc[1])
        itsc1np = np.array(itsc[1])
        if itsc0tp!=(None, None):
            # insert in street1
            name1 = sgmpr[0][6]
            end1 = sgmpr[0][4]
            end2 = sgmpr[0][5]
            end1idx = strts_[name1].index(end1)
            end2idx = strts_[name1].index(end2)
            # index of reference point from which to measure distance
            # the reference point(origin)
            # the subsequent point from origin
            orgidx = min(end1idx, end2idx)
            edidx = max(end1idx, end2idx)
            origin = np.array(strts_[name1][orgidx])
            # construct a list for vertices in a segment
            vtsinseg = []
            for i in range(edidx-orgidx+1):
                vtsinseg.append(np.array(strts_[name1][orgidx+i]))
            # squared distances of vertices
            dtsinseg = []
            for vt in vtsinseg:
                dtsinseg.append(np.sum(np.square(origin-vt)))
            # "squared distance" from origin to the point to be inserted
            dt0 = np.sum(np.square(origin-itsc0np))
            # find the position
            i = 0
            for dt in dtsinseg:
                if dt0<dt:
                    break
                i = i+1
            # check if it is already there
            if itsc0tp not in strts_[name1]:
                # insert in the nearest position
                strts_[name1].insert(orgidx+i, itsc0tp)
            # the second intersection (possible)
            if itsc1tp!=(None, None):
                end1idx = strts_[name1].index(end1)
                end2idx = strts_[name1].index(end2)
                # construct a list for vertices in a segment
                vtsinseg = []
                for i in range(edidx-orgidx+1):
                    vtsinseg.append(np.array(strts_[name1][orgidx+i]))
                # squared distances of vertices
                dtsinseg = []
                for vt in vtsinseg:
                    dtsinseg.append(np.sum(np.square(origin-vt)))
                # "squared distance" from origin to the point to be inserted
                dt1 = np.sum(np.square(origin-itsc1np))
                # find the position
                i = 0
                for dt in dtsinseg:
                    if dt1<dt:
                        break
                    i = i+1
                # check if it is already there
                if itsc1tp not in strts_[name1]:
                    # insert in the nearest position
                    strts_[name1].insert(orgidx+i, itsc1tp)

            # insert in street2
            name2 = sgmpr[1][6]
            end1 = sgmpr[1][4]
            end2 = sgmpr[1][5]
            end1idx = strts_[name2].index(end1)
            end2idx = strts_[name2].index(end2)
            orgidx = min(end1idx, end2idx)
            edidx = max(end1idx, end2idx)
            origin = np.array(strts_[name2][orgidx])
            # construct a list for vertices in a segment
            vtsinseg = []
            for i in range(edidx-orgidx+1):
                vtsinseg.append(np.array(strts_[name2][orgidx+i]))
            # squared distances of vertices
            dtsinseg = []
            for vt in vtsinseg:
                dtsinseg.append(np.sum(np.square(origin-vt)))
            # "squared distance" from origin to the point to be inserted
            dt0 = np.sum(np.square(origin-itsc0np))
            # find the position
            i = 0
            for dt in dtsinseg:
                if dt0<dt:
                    break
                i = i+1
            # check if it is already there
            if itsc0tp not in strts_[name2]:
                # insert in the nearest position
                strts_[name2].insert(orgidx+i, itsc0tp)
            # the second intersection (possible)
            if itsc1tp!=(None, None):
                end1idx = strts_[name2].index(end1)
                end2idx = strts_[name2].index(end2)
                # construct a list for vertices in a segment
                vtsinseg = []
                for i in range(edidx-orgidx+1):
                    vtsinseg.append(np.array(strts_[name2][orgidx+i]))
                # squared distances of vertices
                dtsinseg = []
                for vt in vtsinseg:
                    dtsinseg.append(np.sum(np.square(origin-vt)))
                # "squared distance" from origin to the point to be inserted
                dt1 = np.sum(np.square(origin-itsc1np))
                # find the position
                i = 0
                for dt in dtsinseg:
                    if dt1<dt:
                        break
                    i = i+1
                # check if it is already there
                if itsc1tp not in strts_[name2]:
                    # insert in the nearest position
                    strts_[name2].insert(orgidx+i, itsc1tp)

    return strts_

def gnrtedges(itscs, ends, strts):
    """Generate edges.
    """

    if itscs==[]:
        return None
    else:
        edges = []
        vtxs = itscs+ends
        # go through all the combinations of vertexes
        for i in range(len(vtxs)-1):
            for j in range(i+1, len(vtxs)):
                vt1 = vtxs[i]
                vt2 = vtxs[j]
                # go through all the streets
                for strt in strts.values():
                    # 2 vertexes must be in the same street
                    if vt1 in strt and vt2 in strt:
                        # at least one should be intersection
                        if vt1 in itscs or vt2 in itscs:
                            vt1lb = strt.index(vt1)
                            vt2lb = strt.index(vt2)
                            # 2 vertexes must be near each other
                            if abs(vt1lb-vt2lb)==1:
                                edges.append((i, j))

        # discard the possibly repeated edges
        return list(set(edges))


def main():

    # dictionary for street name and coordinates
    strts = {}

    while True:
        line = sys.stdin.readline()
        if line=='':
            break
        else:
            # regular expression segments
            recmd1 = r'^\s*[ac]' # command1
            recmd2 = r'^\s*r' # command2
            recmd3 = r'^\s*g' # command3
            recmd4 = r'^\s*r\s*$' # command4: remove all
            renm = r'"\s*[a-zA-Z][a-zA-Z\s]*\s*"' # street name
            recd = r'(\(\s*[-]?[0-9]+\.?[0-9]*\s*,\s*[-]?[0-9]+\.?[0-9]*\s*\)\s*){2,}' # coordinates
            # regular expression for line with a or c command
            reln1 = re.compile(recmd1+r'\s+'+renm+r'\s+'+recd+r'$')
            # regular expression for line with r command
            reln2 = re.compile(recmd2+r'\s+'+renm+r'\s*$')
            # regular expression for line with g command
            reln3 = re.compile(recmd3+r'\s*$')
            # regular expression only for 'r'
            reln4 = re.compile(recmd4)

            if reln1.match(line):
                # extract command
                cmd = re.split(r'\s"', line)[0]
                cmd = re.search(r'[ac]', cmd)[0]
                if cmd=='a': # add a street
                    # street name
                    strn = re.search(r'"\s*[a-zA-Z][a-zA-Z\s]*\s*"', line)[0]
                    strn = re.split(r'^"\s*',strn)[-1]
                    strn = re.split(r'\s*"$', strn)[0]
                    # discard repeated spaces
                    if ' ' in strn:
                        words = re.split(r'\s+',strn)
                        strn_ = words[0]
                        for i in range(1, len(words)):
                            strn_ = strn_+' '+words[i]
                        strn = strn_
                    # check if this name exists
                    keys = []
                    for key in strts.keys():
                        keys.append(key.lower())
                    if strn.lower() in keys:
                        sys.stderr.write('Error: {} has already existed.\n'.format(strn))
                        sys.exit(1)
                    else:
                        # street coordinates
                        strcdstr = re.findall(r'[-]?[0-9]+\.?[0-9]*', line) # string list
                        strcds = [] # list of coordinate tuples[(x, y), ...]
                        for i in range(int(len(strcdstr)/2)):
                            strcds.append((float(strcdstr[2*i]), float(strcdstr[2*i+1])))
                        # check repeated end points
                        strcds_ = list(set(strcds))
                        if len(strcds)==len(strcds_):
                            # add name and coordinates in dictionary
                            strts[strn] = strcds
                        else:
                            sys.stderr.write('Error: Should not input repeated points for one street.\n')
                            sys.exit(1)

                else: # cmd=='c' change an existing street
                    # street name
                    strn = re.search(r'"\s*[a-zA-Z][a-zA-Z\s]*\s*"', line)[0]
                    strn = re.split(r'^"\s*',strn)[-1]
                    strn = re.split(r'\s*"$', strn)[0]
                    # discard repeated spaces
                    if ' ' in strn:
                        words = re.split(r'\s+',strn)
                        strn_ = words[0]
                        for i in range(1, len(words)):
                            strn_ = strn_+' '+words[i]
                        strn = strn_
                    # street coordinates
                    strcdstr = re.findall(r'[-]?[0-9]+\.?[0-9]*', line) # string list
                    strcds = [] # list of coordinate tuples[(x, y), ...]
                    for i in range(int(len(strcdstr)/2)):
                        strcds.append((float(strcdstr[2*i]), float(strcdstr[2*i+1])))
                    # specify if the street is existing
                    keys = []
                    for key in strts.keys():
                        keys.append(key.lower())
                    if strn.lower() in keys:
                        # check repeated end points
                        strcds_ = list(set(strcds))
                        if len(strcds)==len(strcds_):
                            # in case of upper and lower cases
                            # appearing at the same time
                            for key in strts.keys():
                                if key.lower()==strn.lower():
                                    del strts[key]
                                    break
                            strts[strn] = strcds
                        else:
                            sys.stderr.write('Error: Should not input repeated points for one street.\n')
                            sys.exit(1)
                    else:
                        sys.stderr.write('Error: "{}" has not been added.\n'.format(strn))
                        sys.exit(1)

            elif reln2.match(line): # cmd=='r' # remove a street
                # street name
                strn = re.search(r'"\s*[a-zA-Z][a-zA-Z\s]*\s*"', line)[0]
                strn = re.split(r'^"\s*',strn)[-1]
                strn = re.split(r'\s*"$', strn)[0]
                # discard repeated spaces
                if ' ' in strn:
                    words = re.split(r'\s+',strn)
                    strn_ = words[0]
                    for i in range(1, len(words)):
                        strn_ = strn_+' '+words[i]
                    strn = strn_
                # specify if the street is existing
                keys = []
                for key in strts.keys():
                    keys.append(key.lower())
                if strn.lower() in keys:
                    for key in strts.keys():
                        if key.lower()==strn.lower():
                            del strts[key]
                            break
                else:
                    sys.stderr.write('Error: "{}" has not been added.\n'.format(strn))
                    sys.exit(1)

            elif reln3.match(line): # cmd=='g'
                if not strts:
                    sys.stderr.write('Error: Cannot generate graph without streets.\n')
                    sys.exit(1)
                elif len(strts)==1:
                    sys.stderr.write('Error: Cannot generate graph with one street.\n')
                    sys.exit(1)
                else:
                    # list storing all the segments of different streets seperately
                    sgmsstrts = []
                    # go through all the streets
                    for nm, cds in strts.items():
                        sgms = []
                        # go through all the coordinates for each street
                        for i in range(len(cds)-1):
                            x1 = cds[i][0]
                            y1 = cds[i][1]
                            x2 = cds[i+1][0]
                            y2 = cds[i+1][1]
                            sgms.append(gnrtseg(x1, y1, x2, y2, nm))
                        sgmsstrts.append(tuple(sgms))
                    sgmsstrts = tuple(sgmsstrts)
                    # dict. {segments pair: intersection, ...}
                    sgpr_itsc = {}
                    # go through all the combinations of different streets
                    for i in range(len(sgmsstrts)-1):
                        for j in range(i+1, len(sgmsstrts)):
                            # go through all the combinations of segments in two streets
                            for k in range(0, len(sgmsstrts[i])):
                                for l in range(0, len(sgmsstrts[j])):
                                    # generate intersections
                                    sgpr_itsc[(sgmsstrts[i][k], sgmsstrts[j][l])] = calc_itsc(
                                        sgmsstrts[i][k], sgmsstrts[j][l])
                    # generate vertexes
                    # (itscs and ends are empty if no intersection,
                    # remove all repeated vertexes in ends)
                    itscs, ends = gnrtvtxs(sgpr_itsc)
                    for itsc in itscs:
                        if itsc in ends:
                            ends.remove(itsc)
                    # insert intersections into streets
                    # (insert nothing if no intersection)
                    strts_ = inrtitscs(sgpr_itsc, strts)
                    # check if there is no intersections
                    # generate edges
                    # (returns None if no intersection)
                    edges = gnrtedges(itscs, ends, strts_)
                    # print vertexes and edges
                    if itscs==[]:
                        print('V = {\n}')
                        print('E = {\n}')
                    else:
                        # output vertex
                        sys.stdout.write('V {}\n'.format(len(itscs+ends)))
                        sys.stdout.flush()
                        # output edges
                        opegs = 'E {'
                        for eg in edges[:-1]:
                            addeg = '<{},{}>,'.format(eg[0], eg[1])
                            opegs = opegs+addeg
                        opegs = opegs+'<{},{}>'.format(edges[-1][0], edges[-1][1])
                        opegs = opegs+'}'
                        sys.stdout.write(opegs+'\n')
                        sys.stdout.flush()

            # clear all stored streets
            elif reln4.match(line):
                strts.clear()

            else:
                sys.stderr.write('Error: bad input.\n')
                sys.exit(1)

    # return exit code 0 on successful termination
    sys.exit(0)


if __name__ == '__main__':
    main()
