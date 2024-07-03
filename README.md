# In-depth Analysis of Densest Subgraph Discovery in a Unified Framework

This repository contains C++ codes and datasets for the paper:

> In-depth Analysis of Densest Subgraph Discovery in a Unified Framework

## Introduction

In this paper, we condcut an in-depth study on sequential densest subgraph discovery (DSD) algorithms. We first propose a unified framework with three modules, namely *graph reduction*, *vertex weight update* ( `VWU`) and *candidate subgraph extract and verify* (`CSV`), which capture the core ideas of all existing algorithms.

Given a graph $G$ and an error threshold $\epsilon$, *graph reduction* aims to locate the DS in a small subgraph;  `VWU` aims to update vertex weights over $T$ iterations; and  `CSV` extracts a candidate subgraph based on vertex weights and verifies if it satisfies the $\epsilon$ error requirement. 

Under this framework, we systematically compare 12 and 7 representative algorithms (including both exact and approximation algorithms) for undirected and directed graphs, respectively.

We conduct comprehensive experiments on both real-world and synthetic datasets and provide an in-depth analysis.

## Environment

The codes of our in-depth study are implemented and tested under the following development environment:

- Hardware : Intel(R) Xeon(R) Gold 6338 CPU @ 2.00GHz and 512GB of memory.
- Operation System : Ubuntu 20.04.5 LTS (GNU/Linux 5.15.0-101-generic x86_64)
## Datasets


We use twelve real datasets from different domains including 6 undirected graphs and 6 directed graphs, which are available on the Stanford Network Analysis Platform, Laboratory of Web Algorithmics, Network Repository, and Konect.


Undirected:


| Dataset         | Category      | $\mid V \mid$   | $\mid E \mid$     | Download link                                                  |
| --------------- | ------------- | -------: | --------: |----------------------------------------------------------------|
| Econ-beacxc (EB) | Economic      | 507     | 42,176    | [Link](https://networkrepository.com/econ-beacxc.php)          |
| DBLP (DP)       | Collaboration | 317,080 | 1,049,866 | [Link](https://snap.stanford.edu/data/com-DBLP.html)           |
| Youtube (YT)    | Multimedia | 3,223,589 | 9,375,374 | [Link](https://snap.stanford.edu/data/com-Youtube.html)        |
|LiveJournal (LJ)|Social|4,036,538|34,681,189| [Link](https://snap.stanford.edu/data/com-LiveJournal.html)    |
|WebBase (WB)|Web|118,142,155|881,868,060| [Link](https://networkrepository.com/web-webbase-2001-all.php) |
|Friendster (FS)|Social|124,836,180|1,806,067,135| [Link](https://snap.stanford.edu/data/com-Friendster.html)     |


Directed:

| Dataset          | Category      |  $\mid V \mid$   | $\mid E \mid$    | Download link                                                |
|------------------| ------------- | -------: | --------: |--------------------------------------------------------------|
| Openflights (OF) |Infrastructure|2,939|30,501| [Link](http://konect.cc/networks/opsahl-openflights/)        |
| Advogato (AD)    |Social|6,541|51,127| [Link](http://konect.cc/networks/advogato/)                  |
| Amazon (AM)      |E-commerce|403,394|3,387,388| [Link](http://konect.cc/networks/amazon0601/)                |
| Bidu-zhishi (BA) |Hyperlink|2,141,300|17,794,839| [Link](http://konect.cc/networks/zhishi-baidu-internallink/) |
| Wiki-en (WE)     |Hyperlink|13,593,032|437,217,424| [Link](http://konect.cc/networks/wikipedia_link_en/)         |
| SK-2005 (SK)     |Web|50,636,154|1,949,412,601| [Link](https://law.di.unimi.it/webdata/sk-2005/)             |



## How to Run the Codes


### A. Code Compilation


After cloning the codes from GitHub, use the following command to compile the codes in the repository :


```sh

cmake CMakeLists.txt

make

```


### B. Command Line Parameters

A general command of our program is like:

```sh

./DensestSubgraph [-option1 value1] [-option2 value2] ...

```

For example, you could run `core-exact` on DP in the following command:

```sh

./DensestSubgraph -path ./data/DP.txt -t u -a e -red k-core -alloc flow-exact -ext flow-exact -ver flow-exact

```

There are a lot of options for you to conduct a thorough evalution among different algorithms:

|Parameters| Value            |Description|
|:---------------|:-----------------|:------------|
|-path| ---              |path to the dataset|
|-t| `u`, `d`         |`u`: undirected, `d`: directed|
|-a| `e`, `a`         |`e`: exact, `a`: approximation|
|-eps| $\epsilon\ge0$   |error threshold for $1+\epsilon$ approximation algorithms|
|-red| refer to B1      |method of *graph reduction*|
|-alloc| refer to B2      |method of `VWU`|
|-ext| refer to B3      |method of *candidate subgraph extraction*|
|-ver| refer to B4      |method of *candidate subgraph verification*|
|-seq| `t`, `f`         |`t`: sequential update strategy, `f`:  simultaneous update strategy|
|-vw| `t`, `f`         |`t`: transform DDS problem into vertex-weighted UDS problem, `f`: do not transform|
|-gamma| $0\le\gamma\le1$ |a parameter that controls the lower bound of binary search|
|-exp| `t`, `f`         |`t`: iteration number grows exponentially, `f`: iteration number is fixed|
|-it| integer, $it\ge1$|fixed iteration number|
|-dc| `t`, `f`         |`t`: apply divide-and-conquer strategy, `f`: do not apply|
|-ra| `t`, `f`         |ablation study on *graph reduction*, `t`: print reduction ratio, `f`: do not print|
|-res| `t`, `f`         |`t`: restrict $xy-core$ in a tight interval, `f`: do not restrict|
|-width| $width\ge1$      |a parameter that controls the tightness of interval|
|-multi| `t`, `f`         |`t`: apply multi-round reduction, `f`: apply single-round reduction|


#### B1. Methods of *Graph Reduction*

|Value|Description|
|--------|--------|
|`k-core`|derive a $k-core$, support UDS algorithms|
|`stable`|derive a stable set|
|`exact-xy-core`|derive an exact $xy-core$, support DDS algorithms|
|`appro-xy-core`|derive an approximate $xy-core$, support DDS algorithms|
|`w-core`|derive an $w^*-core$, support WCoreApp algorithm|


#### B2. Methods of `VWU`

|Value|Description|
|--------|--------|
|`flow-exact`|the `VWU` method of `FlowExact`, `CoreExact`, `DFlowExact`, `DCExact`|
|`fw`|the `VWU` method of `FWExact`, `FWApp`, `DFWExact` and `DFWApp`|
|`fista`|the `VWU` method of `FISTAExact` and `FISTAApp`|
|`mwu`|the `VWU` method of `MWUExact` and `MWUApp`|
|`core-app`|the `VWU` method of `CoreApp`|
|`greedy`|the `VWU` method of `Greedy` and `DGreedy`|
|`greedypp`|the `VWU` method of `Greedy++`|
|`flow-app`|the `VWU` method of `FlowApp`|
|`xy-core-appro`|the `VWU` method of `XYCoreApp`|
|`w-core-appro`|the `VWU` method of `WCoreApp`|


#### B3. Methods of *Candidate Subgraph Extraction* (`CSE`)

|Value|Description|
|--------|--------|
|`flow-exact`|the `CSE` method of `FlowExact`, `CoreExact`, `DFlowExact`, `DCExact`|
|`cp`|the `CSE` method of `FWExact`, `FWApp`,`FISTAExact` ,`FISTAApp`, `MWUExact`, `MWUApp`, `DFWExact` and `DFWApp`|
|`core-app`|the `CSE` method of `XYCoreApp` and `WCoreApp`|
|`greedy`|the `CSE` method of `DGreedy`|


#### B4. Methods of *Candidate Subgraph Verification* (`CSV`)

|Value|Description|
|-------------|--------|
|`flow-exact`|the `CSV` method of `FlowExact`, `CoreExact`, `DFlowExact`, `DCExact`|
|`cp`|the `CSV` method of `FWExact`, `FWApp`,`FISTAExact` ,`FISTAApp`, `MWUExact`, `MWUApp`, `DFWExact` and `DFWApp`|
|`core-app`|the `CSV` method of `CoreApp`|
|`flow-app`|the `CSV` method of `FlowApp`|
|`greedy`|the `CSV` method of `DGreedy`|



### C. Data Download


You can download the datasets from the following Google driven link:


XXXXXXXXXXXXXXXXXX


### D. Experimentation


Our one-click script for reproducibility is comming soon.


[//]: # (### E. Contact)

[//]: # ()
[//]: # ()
[//]: # (If you have any questions about the code or find any errors, please list them in the `issue` or contact us directly by email:)

[//]: # ()
[//]: # ()
[//]: # (`yiyang3@link.cuhk.edu.cn` , `qingshuoguo@link.cuhk.edu.cn` or `yinglizhou@link.cuhk.edu.cn`)
