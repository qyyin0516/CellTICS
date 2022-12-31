import numpy as np
import pandas as pd


def fphi(x):
  zt = (1 - (x + 1).div(x.max(axis=1) + 1, axis=0)).sum(axis=1) / (x.shape[1] - 1)
  zp = (x + 1).div(x.max(axis=1) + 1, axis=0).mul(zt, axis=0)
  return zp


def frho(x):
  zt = (1 - (1 / (x + 1)).mul(x.min(axis=1) + 1, axis=0)).sum(axis=1) / (x.shape[1] - 1)
  zp = (1 / (x + 1)).mul(x.min(axis=1) + 1, axis=0).mul(zt, axis=0)
  return zp


def fmkg(phi, rho, thr1=0.95, thr2=0.9):
  gnm = phi.index.values
  ctpnm = phi.columns.values
  phi = np.array(phi)
  rho = np.array(rho)
  nummkg1 = round((1 - thr1) * phi.shape[0])
  nummkg2 = round((1 - thr2) * phi.shape[0])
  alpha = []
  beta = []
  for i in range(0, phi.shape[1]):
    alpha = np.append(alpha, np.quantile(phi[:, i], thr1))
    beta = np.append(beta, np.quantile(rho[:, i], thr2))
  mkh = []
  mkl = []
  for i in range(0, phi.shape[1]):
    mkh = np.concatenate([mkh, gnm[phi[:, i] >= alpha[i]][0:nummkg1]], axis=0)
    mkl = np.concatenate([mkl, gnm[rho[:, i] >= beta[i]][0:nummkg2]], axis=0)
  mkh = mkh.reshape(nummkg1, phi.shape[1])
  mkl = mkl.reshape(nummkg2, rho.shape[1])
  mkh = pd.DataFrame(mkh, columns=ctpnm)
  mkl = pd.DataFrame(mkl, columns=ctpnm)
  return mkh, mkl


def get_expression(rdata, qdata, rlabel, thrh=0.95, thrl=0.9, normalization=True, marker=True):
    # calculate sum of cell type
    rulabel = rlabel.iloc[:, 0].unique()
    rdt = pd.DataFrame(data=None, columns=None)
    for l in rulabel:
      rdata_l = rdata.iloc[:, rlabel[(rlabel["celltype"] == l)].index.tolist()]
      zs = rdata_l.apply(lambda x: x.sum(), axis=1)
      rdt = pd.concat([rdt, pd.DataFrame(data=zs, columns=[l])], axis=1)

    # normalization
    if normalization:
        rdt_df = rdt
        rdata_df = rdata
        qdata_df = qdata
        rdt = np.array(rdt_df, dtype=np.float32)
        rdata = np.array(rdata_df, dtype=np.float32)
        qdata = np.array(qdata_df, dtype=np.float32)
        rdt = np.divide(rdt, np.sum(rdt, axis=0, keepdims=True)) * 10000
        rdata = np.log2(np.divide(rdata, np.sum(rdata, axis=0, keepdims=True)) * 10000 + 1)
        qdata = np.log2(np.divide(qdata, np.sum(qdata, axis=0, keepdims=True)) * 10000 + 1)

        rdt = pd.DataFrame(data=rdt, columns=rdt_df.columns, index=rdt_df.index)
        rdata = pd.DataFrame(data=rdata, columns=rdata_df.columns, index=rdata_df.index)
        qdata = pd.DataFrame(data=qdata, columns=qdata_df.columns, index=qdata_df.index)

    # match gene ensembl ids between the query data and the reference data
    nm = qdata.index.tolist()
    zg = pd.DataFrame(data=nm, columns=['gene'])
    gid = []
    for i in range(0, len(nm)):
      if nm[i] not in rdata.index.tolist():
        zg.iloc[i] = ''
      else:
        gid = gid + [nm[i]]
    qdata = qdata.loc[gid, :]
    rdata = rdata.loc[gid, :]
    rdt = rdt.loc[gid, :]

    rlabel.index = rdata.columns

    if not marker:
        train_x = rdata
        train_y = rlabel
        test_x = qdata
    else:
        # calculate important genes
        phi = fphi(rdt)
        rho = frho(rdt)
        mkh, mkl = fmkg(phi, rho, thrh, thrl)
        mkg = np.unique(np.append(np.array(mkh).reshape(1, -1), np.array(mkl).reshape(1, -1)))

        # return training set and testing set
        train_x = rdata.loc[mkg, :]
        train_y = rlabel
        test_x = qdata.loc[mkg, :]
    return train_x, test_x, train_y
