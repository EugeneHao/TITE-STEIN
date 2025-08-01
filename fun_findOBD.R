findOBD <- function(tox, eff, pT, qE)
{
  MTD <- sum(tox <= pT)
  if(MTD == 0)
  {
    return(0)
  } else if(max(eff[1:MTD]) < qE)
  {
    return(-1)
  } else
  {
    OBD <- which.max(eff[1:MTD])
    return(OBD)
  }
}





# when u11 = 100, u00 = 0 should provide the same result
findOBD_RDS <- function(tox, eff, pT, qE, u11, u00)
{
  MTD <- sum(tox <= pT)
  if(MTD == 0)
  {
    return(0)
  } else if(max(eff[1:MTD]) < qE)
  {
    return(-1)
  } else
  {
    RDS <- 100 * eff * (1 - tox) + u00 * (1-tox) * (1-eff) + u11 * tox * eff
    OBDlist <- which.max(RDS[1:MTD])
    return(OBDlist)
  }
}
