function [s, state] = state10a(cladeStatus)
    % Return the below tree with or without clade information

    % Root node (1) has age 1000, Adam node is (20)
    %                                             _________________________(3) 3
    %                                           |
    %        __________________________________(2)            _____________(6) 6
    %       |                                   |     ______(5)
    %       |                                   |    |       |_____________(7) 7
    %       |                                   ****(4)
    %       |                                        |           ___(9)* 8
    %       |                                        |_________(8)
    % (20)—(1)                                                  |       __(11) 9
    %       |                                                   |*****(10)
    %       |                                                          |__(12)10
    %       |
    %       |  __________________________(14)* 1
    %       | |
    %       |(13)                                  ________________________(17)4
    %         |                      ____________(16)
    %         |____________________(15)           |________________________(18)5
    %                               |
    %                               |______________________________________(19)2
    % CLADE NAME = c1 ROOTMIN = 450 ROOTMAX = 550 TAXA = 1;
    % CLADE NAME = c2 ROOTMIN = 50  ROOTMAX = 150 TAXA = 8;
    % CLADE NAME = c3 ROOTMIN = 200 ROOTMAX = 500 TAXA = 6 , 7 , 8 , 9 , 10;
    % CLADE NAME = c4                             TAXA = 9 , 10;

    filePath = fullfile('coupling', '+BchooseCoupledMaximal', ...
                        sprintf('state10a-clades%s.mat', cladeStatus));
    state = getfield(load(filePath), 'state');
    s = state.tree;
end
