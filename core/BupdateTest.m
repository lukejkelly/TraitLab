function tests = BupdateTest
    tests = functiontests(localfunctions);
end

function catastropheTest(testCase)
    % mostly checking catastrophe-related outputs
    global BORROWING MCMCCAT NARROW WIDE
    M = 10;
    N = 10;
    for BORROWING = [0, 1]
        for MCMCCAT = [0, 1]
            for mt = [NARROW, WIDE]
                for addClade = [0, 1]
                    for m = 1:M
                        L = 5 + ceil(rand * 5);
                        [state, theta, prior] = dummyState(L, addClade);
                        for n = 1:N
                            newage = [];
                            while isempty(newage)
                                [i, j, k, newage, ~, ncat, cat, loc] = Bchoose(...
                                    state, mt, theta, prior);
                            end
                            [p, q, h] = Bcats.pqh(state, i);
                            % j's parent is now p, h's [sib(i)] parent now pa(pa(i))
                            [nstate, ~, ~] = Bupdate(...
                                state, i, j, k, newage, ncat, cat, loc);

                            assertEqual(testCase, ...
                                        [nstate.tree([j, h, p]).parent], ...
                                        [p, q, k]);
                            if MCMCCAT
                                assertEqual(testCase, nstate.ncat, ncat);
                                assertEqual(testCase, nstate.cat([j, h, p]), ...
                                            structfun(@(x) x, cat));
                                if BORROWING
                                    assertEqual(testCase, ...
                                                {nstate.tree([j, h, p]).catloc}, ...
                                                {sort(loc.j), sort(loc.h), sort(loc.p)});
                                else
                                    assertEqual(testCase, ...
                                                {nstate.tree([j, h, p]).catloc}, ...
                                                {[], [], []});
                                end
                            else
                                assertEqual(testCase, nstate.ncat, state.ncat);
                                assertEqual(testCase, nstate.cat, state.cat);
                                assertEqual(testCase, {nstate.tree.catloc}, ...
                                            {state.tree.catloc});
                            end
                        end
                    end
                end
            end
        end
    end
end

function oldBupdateTest(testCase)
    global BORROWING MCMCCAT NARROW WIDE
    M = 10;
    N = 100;
    MCMCCAT = 0;
    BORROWING = 0;
    for mt = [NARROW, WIDE]
        for addClade = [0, 1]
            for m = 1:M
                L = 10 + ceil(rand * 5);
                [state, theta, prior] = dummyState(L, addClade);
                for n = 1:N
                    newage = [];
                    while isempty(newage)
                        [i, j, k, newage, ~, ncat, cat, loc] = Bchoose(...
                            state, mt, theta, prior);
                    end
                    [nstateOld, UOld, TOld] = oldBupdate(state, i, j, k, newage);
                    [nstateNew, UNew, TNew] = Bupdate(...
                        state, i, j, k, newage, ncat, cat, loc);

                    assertEqual(testCase, nstateOld, nstateNew);
                    assertEqual(testCase, UOld, UNew);
                    assertEqual(testCase, TOld, TNew);
                end
            end
        end
    end
end

function pRootTest(testCase)
    % Bupdate behaviour when p = pa(i) is root
    % Borrowing shouldn't change anything here
    global BORROWING MCMCCAT WIDE
    L = 10;
    M = 50;
    MCMCCAT = 1;
    for BORROWING = [0, 1]
        for m = 1:M
            done = 0;
            % find i who's parent is root
            while ~done
                [state, theta, prior] = dummyState(L, 0);
                bl = getBranchLengths(state);
                for z = find(bl(:)')
                    state.cat(z) = z;
                    if BORROWING
                        state.tree(z).catloc = sort(rand(1, z));
                    end
                end
                r = state.root;
                s = state.tree;
                [i, j, k, newage, ~, ncat, cat, loc] = Bchoose(...
                    state, WIDE, theta, prior);
                if ~isempty(newage) && s(i).parent == r
                    done = 1;
                end
            end
            [p, ~, h] = Bcats.pqh(state, i);
            assertEqual(testCase, p, r);

            jc = state.cat(j);
            hc = state.cat(h);
            pc = state.cat(p);
            assertEqual(testCase, jc, j);
            assertEqual(testCase, hc, h);
            assertEqual(testCase, pc, 0);

            % h is the root in new state
            assertEqual(testCase, cat.j + cat.p, jc);
            assertEqual(testCase, cat.h, 0);

            [nstate, ~, ~] = Bupdate(state, i, j, k, newage, ncat, cat, loc);
            assertEqual(testCase, nstate.cat([j, h, p])', [cat.j, cat.h, cat.p]);
            if BORROWING
                assertLength(testCase, nstate.tree(j).catloc, cat.j);
                assertLength(testCase, nstate.tree(h).catloc, cat.h);
                assertLength(testCase, nstate.tree(p).catloc, cat.p);
            end
            inds = setdiff(find(bl(:)'), [j, h, p]);
            assertEqual(testCase, nstate.cat(inds), state.cat(inds));
        end
    end
end

function jRootTest(testCase)
    % Bupdate behaviour when j is root
    % Borrowing shouldn't change anything here
    global BORROWING MCMCCAT WIDE
    L = 10;
    M = 100;
    MCMCCAT = 1;
    for BORROWING = [0, 1]
        for m = 1:M
            done = 0;
            % find i who's parent is root
            while ~done
                [state, theta, prior] = dummyState(L, 0);
                bl = getBranchLengths(state);
                for z = find(bl(:)')
                    state.cat(z) = z;
                    if BORROWING
                        state.tree(z).catloc = sort(rand(1, z));
                    end
                end
                r = state.root;
                [i, j, k, newage, ~, ncat, cat, loc] = Bchoose(...
                    state, WIDE, theta, prior);
                if ~isempty(newage) && j == r
                    done = 1;
                end
            end
            [p, ~, h] = Bcats.pqh(state, i);

            jc = state.cat(j);
            hc = state.cat(h);
            pc = state.cat(p);
            assertEqual(testCase, jc, 0);
            assertEqual(testCase, hc, h);
            assertEqual(testCase, pc, p);

            % p is the root in new state
            assertGreaterThanOrEqual(testCase, cat.j, 0);
            assertEqual(testCase, cat.h, hc + pc);
            assertEqual(testCase, cat.p, 0);

            [nstate, ~, ~] = Bupdate(state, i, j, k, newage, ncat, cat, loc);
            assertEqual(testCase, nstate.cat([j, h, p])', [cat.j, cat.h, cat.p]);
            if BORROWING
                assertLength(testCase, nstate.tree(j).catloc, cat.j);
                assertLength(testCase, nstate.tree(h).catloc, cat.h);
                assertLength(testCase, nstate.tree(p).catloc, cat.p);
            end
            inds = setdiff(find(bl(:)'), [j, h, p]);
            assertEqual(testCase, nstate.cat(inds), state.cat(inds));
        end
    end
end

function [state, theta, prior] = dummyState(L, addClade)
    global ROOT MCMCCAT
    theta = rand * 1e-2;
    s = ExpTree(L, theta);
    if addClade
        clade = synthclades(s, 1, 2, 1 - rand^3);
        rootmax = (1 + rand) * s([s.type] == ROOT).time;
        prior = unitTests.clade2prior(clade, rootmax);

        state = unitTests.dummyState(s);
        state = UpdateClades(state, [state.leaves, state.nodes], ...
                             size(prior.clade, 2));
    else
        state = unitTests.dummyState(s);
        prior = struct('isclade', 0);
    end
    if MCMCCAT
        for c = 1:poissrnd(4)
            state = AddCat(state);
        end
    end
    state.rho = state.ncat / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    global VARYRHO
    unitTests.setupOnce(testCase);
    testCase.TestData.VARYRHO = VARYRHO;
    VARYRHO = 1;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

% Old Bupdate for comparison
function [state,U,TOPOLOGY] = oldBupdate(state,i,j,k,newage)
    % Catastrophe behaviour has changed so we only expect this to match the new
    % function when BORROWING is false and no catastrophes
    % Old comments deleted for brevity
    global OTHER ANST ROOT BORROWING

    s=state.tree;
    Root=state.root;

    iP=s(i).parent;
    PiP=s(iP).parent;
    CiP=s(iP).child(OTHER(s(i).sibling));

    s(PiP).child(s(iP).sibling)=CiP;
    s(CiP).parent=PiP;
    s(CiP).sibling=s(iP).sibling;

    s(j).parent=iP;
    s(k).child(s(j).sibling)=iP;
    s(iP).sibling=s(j).sibling;
    s(iP).child(OTHER(s(i).sibling))=j;
    s(j).sibling=OTHER(s(i).sibling);
    s(iP).parent=k;

    if BORROWING
        if iP == Root
            s(Root).type = ANST;
            Root = CiP;
            s(CiP).type = ROOT;

            state.cat(iP) = state.cat(CiP);
            s(iP).catloc = s(CiP).catloc;

            state.cat(CiP) = 0;
            s(CiP).catloc = [];
        elseif j == Root
            s(Root).type = ANST;
            Root = iP;
            s(iP).type = ROOT;

            state.cat(j) = state.cat(iP);
            s(j).catloc = s(iP).catloc;

            state.cat(iP) = 0;
            s(iP).catloc = [];
        end
    else
        if iP==Root
            s(Root).type=ANST;
            Root=CiP;
            s(Root).type=ROOT;
            state.ncat=state.ncat-state.cat(Root);
            state.cat(Root)=0;
        elseif j==Root
            s(Root).type=ANST;
            Root=iP;
            s(Root).type=ROOT;
            state.ncat=state.ncat-state.cat(Root);
            state.cat(Root)=0;
        end
        state.cat(CiP)=state.cat(CiP)+state.cat(iP);
        state.cat(iP)=0;

        state.cat(iP)=binornd(state.cat(j),(newage-s(j).time)/(s(k).time-s(j).time));
        state.cat(j)=state.cat(j)-state.cat(iP);
    end

    oldage = s(iP).time;
    s(iP).time = newage;

    U=above([PiP,iP],s,Root);

    state.tree=s;
    state.root=Root;

    if any([iP, PiP, CiP, j]==Root)
        state.length=TreeLength(state.tree,Root);
    else
        state.length=state.length-oldage+newage;
    end
    TOPOLOGY=1;
end
