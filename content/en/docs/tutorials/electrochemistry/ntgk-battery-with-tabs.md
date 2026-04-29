---
title: "NTGK Battery With Tabs"
linkTitle: "NTGK Battery With Tabs"
weight: 1
description: >
  Thermal simulation of a battery cell using the NTGK electrochemical model.
---

**Repository path:** `tutorials/electrochemistry/NTGKBatteryWithTabs`

## Overview

This case demonstrates the coupled electrochemical–thermal simulation of a
lithium-ion battery cell including current collector tabs. The electrochemistry
is described by the Newman–Tiedemann–Gu–Kim (NTGK) model, coupled to a heat
transport equation across the cell regions.

## Running the case

```bash
cd $FOAM_RUN/../multiPhysicsFoam/tutorials/electrochemistry/NTGKBatteryWithTabs
./Allrun
```
